#!/usr/bin/env python3

"""
#-------------------------------------
#Program Name: volvox_scan
#Version: Release Version 1.0
#Date Written: 08/15/2025 
#E-mail Contact: volvoxscan@gmail.com
#-------------------------------------

#Program Description:

The volvox_scan DNA mutation clustering software is a standalone program with code written entirely using Python language (version 3). The program predicts likely false-negative DNA mutations and likely related mutants from a sequenced population of independently mutagenized individuals. The program utilizes variant call format (VCF) files containing unfiltered, initial variant-calls from single-individual genotyping calls. The program does not currently support VCF files for joint-genotyping calls. By default the program requires no input parameters and searches for any collection or mitxure of commpressed (zipped) and/or uncompressed VCF files in the current folder (or directory) where the program is executed or invoked. On successful program analysis completion the results will be located in a "Results" folder. The program is platform-independent and has been tested successfully on Windows, UNIX/Linux and Mac OS.  

#Example Runs:
#python3 volvox_scan.py --snps -d ../data/both/dupl 
#python3 volvox_scan.py --all -q 10 -d ../data/both/dupl
#python3 volvox_scan.py --all -q 10 -d ../data/v2.1_255/ -n 5 -m 250 
#python3 volvox_scan.py --all -q 10 --uniq -d ../data/v2.1_255/ -n 5 -m 250 >& xn.txt

# Dictionary Descriptions:
#
# -----------------------------------------------------------------------
#| Dictionary |                 Key                 |       Value        |
#|-----------------------------------------------------------------------|
#| allvar     | SampleID + Chrom (or Scaffold) Name | Entire VCF Line    |
#|            | + Genome Co-ordinate                |                    |
# -----------------------------------------------------------------------
#| posns      | Chrom (or Scaffold) Name            | List of Sample IDs |
#|            | + Genome Co-ordinate                |                    |
#|-----------------------------------------------------------------------|
#| counts     | Chrom (or Scaffold) Name            |  Sample Count      |
#|            | + Genome Co-ordinate                |                    |
# -----------------------------------------------------------------------
#| totals     | Concatenated List of Sample IDs     |  Number of Shared  |
#|            | Separated By Pipes ("|")            |  Mutations         |
# -----------------------------------------------------------------------
"""

try:
    FILES="volvox_scan_vcf_files.txt"

    import os, gzip, argparse, shutil, platform, datetime

except ImportError as e:
    print("ERROR: A Required Python Module is Not Installed.")
    print("ERROR:",e)
    exit()


def main():
    DATE        = volvox_scan_printDate()   #Record date and time at start of analysis 
    MAXCLUSTPOP = 5         #Default maximum number of samples in a retained cluster
    MINCLUST    = 250       #Default minimum number of shared mutations in a cluster
    MINQUAL     = 10        #Default minimum variant quality score for retained mutations
    PATHS       = os.getcwd()
    PREFIX      = ""        #Default Prefix for Analysis Summary Filenames
    allvar      = {}        #Dictionary of All Variants detected in all VCF Files
    posns       = {}        #Dictionary of Distinct Mutation Positions and Corresponding list of Sample IDs
    counts      = {}        #Dictionary of Distinct Mutation Positions and Corresponding Sample Counts
    totals      = {}        #Dictioanry of Sample ID List and Total Number of Shared Mutations
    cluster     = []
    overFlow    = []        # List of additional (unused) multiple variants found at the sample genome positions in a single sample

    #Detect Operating System:
    OSTYPE = platform.system()
    #print("Detected OS:",OSTYPE)

    #Detected Linux or Mac Operating System:
    if OSTYPE != "Windows":

        parser = argparse.ArgumentParser(description="Detection of Clusters of Likely Related Mutants in a Sequenced Population")
        parser.add_argument('-d','--dir',help='Path For Directory Containing VCF Files', type=str,default=PATHS)
        parser.add_argument('-p','--prefix',help='Prefix For Output File Names',default=PREFIX,type=str)
        parser.add_argument('-q','--qual',help='Minimum Variant Quality for Mutation Spectrum, Heterozygosity and Common Variant Filtering Analyses', type=float, default=MINQUAL)
        parser.add_argument('-m','--min',help='Minimum No. of Shared Variants to Report a Cluster', type=int, default=MINCLUST)
        parser.add_argument('-n','--maxpop',help='Maximum No. of Individuals in a Cluster', type=int, default=MAXCLUSTPOP)
        parser.add_argument("--snps",help='Use SNPs Only in Clustering Analysis',action="store_true")
        parser.add_argument("--all",help='Use All Variants in Clustering Analysis' ,action="store_true")
        parser.add_argument("--uniq",help='Perform Common Variants Subtraction',action="store_true")

        #Retrieve Parameters from Command Line Arguments:
        args        = parser.parse_args()
        DIR         = args.dir
        PREFIX      = args.prefix
        MINCLUST    = args.min
        MINQUAL     = args.qual
        MAXCLUSTPOP = args.maxpop

        if args.all and args.snps:
            print("\n\nERROR: Choose Either --snps or --all Option and NOT BOTH\n\n")
            exit()

        else:
            if args.all == True:        #Use All Variants in Clustering Analysis
                SNPS = 0
            elif args.snps == True:     #Use SNPs Only    in Clustering Analysis
                SNPS = 1
            else:
                SNPS = 1

        if args.uniq == True:           #Perform Common Variants Subtraction in Analysis
                UNIQS = 1
        else:
                UNIQS = 0

    #Detected Windows Operating System:
    else:   #Windows OS:

        #Prompt User for Parameter Values:
        promptStr   = volvox_scan_getPromptStrings()
        DIR         = input("Enter File Path to VCFs Directory: ")
        PREFIX      = input("Enter Prefix String For Output Files: ")
        MINCLUST    = int(volvox_scan_getParam(MINCLUST, promptStr[0]))
        MINQUAL     = float(volvox_scan_getParam(MINQUAL, promptStr[1]))
        MAXCLUSTPOP = int(volvox_scan_getParam(MAXCLUSTPOP, promptStr[2]))
        SNPS        = volvox_scan_getValues(1,promptStr[3])
        UNIQS       = volvox_scan_getValues(0,promptStr[4])

        if len(DIR.strip()) == 0:
            DIR = os.getcwd()

    #END OF ELSE WINDOWS:

    SEQS="UNKNOWN"
    if SNPS == 1:
        SEQS="SNPS ONLY"
    else:
        SEQS="SNPS And INDELS"
    
    DIR  = DIR.strip()
    FULL = os.path.expanduser(DIR)   #Instances of tilde in path
    FULL = os.path.abspath(FULL)     #Generate Full Directory Path
   
    #Display All Parameters Selected:
    print("Detected OS    : ", OSTYPE)
    print("DIR            : ", FULL)
    print("PREFIX         : ", PREFIX)
    print("MINCLUST       : ", MINCLUST)
    print("MINQUAL        : ", MINQUAL)
    print("MAXCLUSTPOP    : ", MAXCLUSTPOP)
    print("VARIATION TYPE : ", SEQS)
    #print("\n\n")

    try:
        CURRDIR = os.getcwd()
        RESULTS = os.path.join(CURRDIR,"Results")

        #Check and Caution User if a (previous) "Results" folder already Exists:
        if (os.path.exists(RESULTS) and os.path.isdir(RESULTS)):

            response = input("An Existing Results Folder Exists. Delete Folder (Y/N)?: ")

            if response.upper() == "Y":
                print("Cleaning Up Results Folder ....")
                shutil.rmtree(RESULTS)
            else:
                print("Aborting Program Execution ...")
                exit()

        #Create Folders For Analysis Reporting:
        os.mkdir("Results")
        RESULTS = os.path.abspath("Results")
        os.chdir(RESULTS)
    
        os.mkdir("Summary")
        SUMMARY = os.path.abspath("Summary")

        os.mkdir("Shared_VCFs")
        SHARED = os.path.abspath("Shared_VCFs")

        os.mkdir("Clusters")
        CLUSTERS = os.path.abspath("Clusters")

        if UNIQS == 1:
            os.mkdir("Uniques_VCFs")
            UNIQUES = os.path.abspath("Uniques_VCFs")

    except OSError as err:
        print("Directory Creation Error: ", err)
        exit()

    uniqCnt = -999       #Initialize Unique Count
    dataStat=[0,0,0]     #Read VCFs Statistics

    #Retrieve and Read VCF Files in Source Directory and Create Analysis Data Structures
    num_samples = volvox_scan_readVCFs(FULL,SNPS,allvar,posns,counts,overFlow,dataStat)
    
    if num_samples < 2:
        print("Population Size is Invalid: ", num_samples)
        exit()


    #Generate Report on Distinct Mutation Positions and their Corresponding List of Origin Samples IDs 
    volvox_scan_create_genome_position_stats(num_samples,posns,counts,PREFIX)

    #Predict Clusters of Samples with Large Number of Shared Mutations:
    mutStat = volvox_scan_predict_clusters(num_samples,allvar,posns,counts,MINCLUST,MAXCLUSTPOP,PREFIX)

    #Report Mutation Spectrum and Computer Heterozygosity for likely False-Negative Shared Mutations:
    volvox_scan_spectrum_zygosity(MINQUAL,MAXCLUSTPOP,PREFIX)

    #Sort the Variants in Shared Variants VCF Files:
    volvox_scan_sortVCF()


    #Perform Common Variants Subtraction:
    if UNIQS == 1:
        print("Filtering Common Variants ...")
        uniqCnt = volvox_scan_filterVCFs(FILES,SNPS,counts,MINQUAL)
        print("Total Number of Filtered Variants or Uniques: ",uniqCnt)
   
    #Report Skipped Variants from Multiple Variants Detected at Same Position in a Single Sample:
    volvox_scan_printOverFlow(overFlow,PREFIX)

    #Generate Summary Clustering Analysis Report:
    volvox_scan_printStat(dataStat,mutStat,uniqCnt,overFlow,MINQUAL,PREFIX,DATE)

    Files = os.listdir(RESULTS)

    #Organize Analysis Files into Respective Sub-Directories:
    for fname in Files:
        SOURCE = os.path.join(RESULTS,fname)

        if (os.path.exists(SOURCE) and os.path.isfile(SOURCE)):
            if fname.endswith("-shared.vcf"):
                DEST = os.path.join(SHARED,fname)
            elif fname.endswith("-uniq.vcf") or fname.endswith("-uniq.vcf.gz"):
                DEST = os.path.join(UNIQUES,fname)
            elif "cluster-" in fname and fname.endswith(".txt"):
                DEST = os.path.join(CLUSTERS,fname)
            else:
                DEST = os.path.join(SUMMARY,fname)

            os.rename(SOURCE,DEST)
        
    return  #End of Main

def volvox_scan_readVCFs(DIR,SNPS,allvar,posns,counts,overFlow,dataStat):
    #Detect and Read all VCF files in the input source Directory provided by user
    #Read variants in all detected VCF file and create lookup dictionaries
    #Also generates counts for all variant categories 
    #The default directory is the current user folder

    if not (os.path.exists(DIR) and os.path.isdir(DIR)):
        print("ERROR: Directory Not Found: ", DIR)
        exit()

    fileCnt = 0
    byteCnt = 0
    #Retrieve list of files and subfolders in source directory:
    Files   = os.listdir(DIR) 

    try:
        fp = open(FILES,"w")
        #print(FILES)

    except (OSError, IOError) as e:
        print("ERROR: File Access Error!!!")
        print("ERROR:",e)
        print("Could NOT create File %s: " % (FILES))
        exit()

    for fname in Files:
        pname = DIR+os.sep+fname

        #Found a Zipped VCF File:
        if pname.endswith(".vcf.gz") and os.path.isfile(pname):
            
            fp.write(pname+"\n") 
            #Retrieve Sample ID Name From File Path:
            sample = volvox_scan_getSampleID(pname)
            print("Processing Sample: ", sample)
            fileCnt += 1
            byteCnt += os.path.getsize(pname)
            print("Opening zipped File: ", pname, os.path.getsize(pname),"bytes")
            with gzip.open(pname,"rt") as fpz:
   
                #Read VCF File and update the variants dictionary data structures:
                dStat = volvox_scan_databases(fpz,sample,SNPS,allvar,posns,counts,overFlow)

                dataStat[0] += dStat[0]     #Update Total Variants Count
                dataStat[1] += dStat[1]     #Update Total SNPs Count
                dataStat[2] += dStat[2]     #Update Total INDELs Count


            fpz.close()


        #Found an Unzipped VCF File:
        elif pname.endswith(".vcf") and os.path.isfile(pname):
    
            fp.write(pname+"\n") 
            #Retrieve Sample ID Name From File Path:
            sample = volvox_scan_getSampleID(pname)
            print("Processing Sample: ", sample)
            fileCnt += 1
            byteCnt += os.path.getsize(pname)
            print("Opening unzipped File: ", pname, os.path.getsize(pname),"bytes")
            with open(pname,"r") as fpz:

                #Read VCF File and update the variants dictionary data structures:
                dStat = volvox_scan_databases(fpz,sample,SNPS,allvar,posns,counts,overFlow)

                dataStat[0] += dStat[0]     #Update Total Variants Count
                dataStat[1] += dStat[1]     #Update Total SNPs Count
                dataStat[2] += dStat[2]     #Update Total INDELs Count


            fpz.close()

        else:   #Other Files and Sub-Folders in Source Directory: 
            print("Sneaky File or SubFolder: ", pname)
            pass


    dataStat.append(fileCnt)    #Insert Total Number of Files
    dataStat.append(byteCnt)    #Insert Total Number of Bytes in Files
    num_samples = fileCnt


    fp.close()

    return num_samples  #End of volvox_scan_readVCFs

def volvox_scan_getSampleID(pname):
    #Extract Sample ID Name from VCF file path:

    xname  = pname
    if pname.endswith(".vcf.gz"):
        xname  = xname.split(".vcf.gz")

    elif pname.endswith(".vcf"):
        xname  = xname.split(".vcf")
    else:
        print("ERROR: Unknown VCF Filename Extension: ", pname)
        exit()

    xname  = xname[0].split(os.sep)
    num    = len(xname)
    num   -= 1
    sample = xname[num]

    return sample   #End of volvox_scan_getSampleID

def volvox_scan_getSampleName(pname):
    #Extract Sample ID Name from analysis filename:

    sample = ""
    xname  = pname
    xname  = xname.replace("-shared.vcf","")
    xname  = xname.replace("cluster-","")

    #Search for second underscore if filename:
    k      = 0
    LEN    = len(xname)

    for i in range(LEN):
        if xname[i] == '_':
            k += 1

        if k == 2:
            break

    i += 1      # Get Next Position After "_"      
    if i >= LEN:
        print("\n\nERROR: File Name Index Out of Range:",pname,"\n\n")
        exit()
    else:
        sample = xname[i:]

    return sample   #End of volvox_scan_getSampleName

def volvox_scan_databases(fp,sample,SNPS,allvar,posns,counts,overFlow):
    #Read variants information in VCF file and updates the variant counts and dictionaries

    varStat   = 0
    snpStat   = 0
    indelStat = 0

    while True:

        line = fp.readline()

        if line == "":
            break

        if line[0] =="#":
            continue
    
        varStat   += 1   #Increment Variants count

        line   = line.strip()
        word   = line.split("\t")
        sample = sample.strip()
        posn   = word[0].strip()+":"+word[1].strip()
        ref    = word[3].strip()
        allele = word[4].strip()

        if "INDEL" in line:
            indelStat += 1   #Increment Indel count 
        
        elif (len(ref) == 1 and len(allele) == 1):
            snpStat   += 1   #Increment SNP count
   
        if (SNPS == 1) and  (len(ref) != 1 or len(allele) != 1):
            continue    #skip

        #Create Dictionary Key:
        seq    = sample+":"+posn   #Prepend sample ID to genome position
        
        sline  = ""
        #Create Dictionary Value:
        for i in range(0,len(word)):
            sline += (word[i]+"\t")
        sline      = sline.strip()

        #Insert Dictionary Record:
        #Check for Multiple Variants at Same Position in a single Individual

        result = allvar.get(seq, "None")

        if result == "None":
            allvar[seq] = sline
        else:
            over = sample+"\t"+sline
            overFlow.append(over)   #Multiple variants at same genomic position

        #Check Position dictionary:
        result = posns.get(posn, "None")

        if result == "None":                # New Genomic Position
            posns[posn]   = sample          # Create Position Entry
            counts[posn]  = 1               # Set Distinct Count to 1

        else:                               # Shared Position Exists

            if sample+"|" in result+"|":
                pass

            else:                               # New Distinct Sample ID
                posns[posn]  += ("|"+sample)    # Append Sample ID to growing list
                counts[posn] += 1               # Update Shared Position Sample Count

    dStat = [varStat,snpStat,indelStat]

    return  dStat   #End of volvox_scan_databases

def volvox_scan_create_genome_position_stats(num_samples, posns, counts,prefix):

    if len(posns) != len(counts) :
        print("ERROR: Dictionary Inconsistencies!!!")
        exit()
    else:
        tally = {}      # Tally Dictionary

        # Initialize Tally Table:
        for i in range(1,num_samples+1):
            tally[i] = 0

        pname = "Position_Cnts.txt"
        sname = "Position_Distribution.txt"

        if len(prefix) > 0:
            pname = prefix.strip()+"-" + pname
            sname = prefix.strip()+"-" + sname

        try:
            fp = open(pname,"w")

        except (OSError, IOError) as e:
            print("ERROR: File Access Error!!!")
            print("ERROR:",e)
            print("Could NOT create File %s: " % (pname))
            exit()

        try:
            fp1 = open(sname,"w")

        except (OSError, IOError) as e:
            print("ERROR: File Access Error!!!")
            print("ERROR:",e)
            print("Could NOT create File %s: " % (sname))
            exit()

        fp.write("Genome Position\tOccurrences\tSample IDs\n")
        for posn in posns:
            result = counts.get(posn, "None")

            if result != "None":
                fp.write("%s\t%d\t%s\n" % (posn,result,posns[posn]))

                if result <= num_samples:
                    tally[result] += 1
                else:
                    pass    #Verify Python pass

            else:
                print("ERROR: Position NOT Found!!!")
                exit()

        fp1.write("Occurrences\tNumber of Distinct Variant Positions\tTotal Number of Variants\n")
        for i in range(1,num_samples+1):
                fp1.write("%d\t%d\t%d\n" % (i,tally[i], i*tally[i]))

        fp.close()
        fp1.close()

    return  #End of volvox_scan_create_genome_position_stats

def volvox_scan_predict_clusters(num_samples, allvar, posns, counts,minCnt,maxClust,prefix):

    new_counts = {}     #Dictionary of Non-Zero and Non-Singleton Position Counts
    cname = "Cluster_Sizes.txt"

    if len(prefix) > 0:
        cname = prefix.strip()+"-" + cname

    try:
        fp = open(cname,"w")
        fp.write("SAMPLE IDs\tCLUSTER ID\tSAMPLE COUNT\tSHARED MUTATIONS\n")

    except (OSError, IOError) as e:
        print("ERROR: File Access Error!!!")
        print("ERROR:",e)
        print("Could NOT create File %s: " % (cname))
        exit()
   
    #Create New Dictionary of Positions with Minimum of two Sample IDs
    for i in counts:

        if counts[i] > 1:
            new_counts[i] = counts[i]   #Get Entries with minimum of Two Sample IDs


    clustCnt     = 1        #Final Cluster ID Counts for clusters with minimum no. of shared SNPs.
    clustCounter = 1        #All Clusters ID Counts
    loc          = []       #List of Genomic Positions shared by clusters of size N
    cluster      = []       #List of Candidate Cluster Sample IDs containing N samples 
    new_cluster  = []       #List of Candidate Cluster Sample IDs containing N samples 
    totals       = {}       #Dictionary of Sample ID List and Number of Shared Mutations

    for i in new_counts:

        # Get all Genomic Positions shared by i-samples:
        # Create a List of sample IDs containing i-samples
        # Create a List of Positions shared by the i-samples

        tword = posns[i].split("|")     #Retrieve Sample ID Names
        tword.sort()                    #Sort List of Sample IDs
        sword = ""
        for n in range(len(tword)):     #Re-Make Sorted Sample ID List
            sword +=(tword[n] + "|")
        tup   = (i,sword)               #Create tuple of Positions and Sample ID List
        new_cluster.append(tup)

    #Sort the List of tuples using the sample list field:

    new_cluster.sort(key=lambda x:x[1]) #Sort List of tuples by Using Sample ID List as Sort Key

    new_count = len(new_cluster)
    match     = []      #List of Concatenated (Using tabs i.e. "\t") Positions and Sample ID List
    sites     = []      #List of Concatenated (Using colons i.e. ":") Positions and Sampel ID List
    n         = 0       #Counter Variable for List of Tuples
    cnt       = 0       #Shared Mutations Counter
    prev      = ""      #Current List of Sample IDs
    
    passSizes     = 0    #Sum of Retained Clusters Sizes 
    passShared    = 0    #Sum of Retained Distinct Shared Mutations
    passMutations = 0    #Sum of Retained All Shared Mutations (i.e. passSizes * passShared)
   
    #Create Clusters:
    while n < new_count:

        if prev != new_cluster[n][1]:   #Compare Current and New Row of Sample ID Lists
                
            if cnt > 0:
                totals[sample]= cnt     #Insert Sample List and Count of Shared Mutations
                
                #Retained Cluster must contain:
                #At least minCnt Number of Shared Mutations
                #And at most maxClust Number of Samples or Mutants
                if cnt >= minCnt and size <= maxClust:  #Check Criteria for Retained Cluster
                    clustCnt += 1       #Update Cluster ID Counter

                    #Generate a Shared Mutations VCF File for each Sample in Cluster:
                    volvox_scan_cluster_write(allvar,match,sites, clustCnt, size, cnt)
            
                    passSizes     += size           #Update Sum of Retained Clusters Sizes 
                    passShared    += cnt            #Update Sum of Retained Distinct Shared Mutations
                    passMutations += (size * cnt)   #Update Sum of Retained All Shared Mutations
   
            #Retrieve Details for New List of Sample IDs
            prev   = new_cluster[n][1]
            end    = len(prev) - 1      #Exclude last pipe character
            sample = prev[:end]         #Retrieve Sample ID List For Next Candidate Cluster
            word   = sample.split("|")
            size   = len(word)
            match  = []
            sites  = []
            cnt    = 0

        #Create List of Search Keys for retrieving data for Later Creation of Shared VCF Files:
        for j in word:
            hit = j+"\t"+new_cluster[n][0]  #Create Sample and Position Record
            match.append(hit)
            hit = j+":" +new_cluster[n][0]  #Create Search Key for Querying allvar Dictionary
            sites.append(hit)
                            
        cnt += 1    #Update Shared Mutations Counter

        n += 1      #Update Index to next record in List of tuples

    #Handle Last Cluster After Termnating Clustering While-Loop Above:
    if cnt > 0:
        totals[sample]= cnt
                
        if cnt >= minCnt and size <= maxClust:
            clustCnt += 1
            volvox_scan_cluster_write(allvar,match,sites, clustCnt, size, cnt)
            passSizes     += size           #Update Sum of Retained Clusters Sizes 
            passShared    += cnt            #Update Sum of Retained Distinct Shared Mutations
            passMutations += (size * cnt)   #Update Sum of Retained All Shared Mutations
    
            
    
    SumSizes     = 0    #Sum of All Clusters Sizes  (All = Both Retained and Non-Retained Clusters)
    SumShared    = 0    #Sum of All Distinct Shared Mutations
    SumMutations = 0    #Sum of All Redundant Shared Mutations
   
    #Report or Create a File of All Clusters Found
    for i in totals:
        iword = i.split("|")
        size  = len(iword)

        #Note: Prefix "clust-" is used for All Cluster IDs
        #Include the number of samples in Cluster (size) and number of shared mutations:
        fp.write("%s\tclust-%d_%d\t%d\t%d\n" % (i, clustCounter, size, size, totals[i]))
        clustCounter +=1    #Cluster ID Counter

        SumSizes     += size                #Update Sum of Clusters Sizes
        SumShared    += totals[i]           #Update Sum of Distinct Shared Mutations
        SumMutations += (size * totals[i])  #Update Sum of All Shared Mutations
    
    fp.close()

    mutStat = [clustCounter-1, SumSizes, SumShared, SumMutations, clustCnt-1, passSizes, passShared, passMutations]

    return mutStat  #End of volvox_scan_predict_clusters

def volvox_scan_cluster_write(allvar,match,sites, clustCnt, size, cnt):

    #Note: Prefix "cluster-" is used for Retained Cluster IDs
    #Include the number of samples in Cluster (size) and number of shared mutations in filename:
    Clust  = "cluster-"+str(clustCnt)+"_"+str(size)+"_"
    filename  = Clust+str(cnt)+".txt"
    clustCnt += 1

    try:
        fp1 = open(filename,"w")

    except (OSError, IOError) as e:
        print("ERROR: File Access Error!!!")
        print("ERROR:",e)
        print("Could NOT create File %s: " % (filename))
        exit()

    #For each Sample in Cluster Report all Sample ID and Shared Mutation Positions
    for name in match:
        fp1.write("%s\n" %(name))

    fp1.close()

    #For each Sample in Cluster Report all VCF Info for al Shared Mutations:
    for hit in sites:

        word = hit.split(":")
        j    = word[0].strip()

        result = allvar.get(hit, "None")    #Search for VCF Line Info corresponding to Position in Sample

        if result != "None":
            filename = Clust+j+"-shared.vcf"
                            
            try:
                fps = open(filename,"a")

                if os.stat(filename).st_size == 0:
                    #CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  Sample
                    fps.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n")

            except (OSError, IOError) as e:
                print("ERROR: File Access Error!!!")
                print("ERROR:",e)
                print("Could NOT create File %s: " % (filename))
                exit()

            fps.write("%s\n" % (result))    #Write VCF File Line into Shared VCF File for a Sample
            fps.close()
        else:
            print("Error: Sample and Genomic Location Not Found", hit)
            exit()

    return      #End of volvox_scan_cluster_write

def volvox_scan_spectrum_zygosity(QUAL,MAXCLUSTPOP,prefix):

    outname = "Mutation_Spectrum.txt"
    zygname = "Percent_Heterozygosity.txt"
    zygsort = "Sorted_Zygosity.txt"

    if len(prefix) > 0:
        outname = prefix.strip()+"-" + outname
        zygname = prefix.strip()+"-" + zygname
        zygsort = prefix.strip()+"-" + zygsort

    try:
        fp1 = open(outname, "w")
        fp1.write("Shared Sample ID\tG:C->A:T(%)\tG:C->T:A(%)\tG:C->C:G(%)\tA:T->G:C(%)\tA:T->C:G(%)\tA:T->T:A(%)\tOthers(%)\t")
        fp1.write("G:C->A:T\tG:C->T:A\tG:C->C:G\tA:T->G:C\tA:T->C:G\tA:T->T:A\tOthers\tTotal\n")

    except (OSError, IOError) as e:
        print("ERROR: File Access Error!!!")
        print("ERROR:",e)
        print("Could NOT create File %s: " % (outname))
        exit()

    try:
        fp2 = open(zygname, "w")
        fp2.write("Shared Sample ID\tHeterozygosity(%)\n")

    except (OSError, IOError) as e:
        print("ERROR: File Access Error!!!")
        print("ERROR:",e)
        print("Could NOT create File %s: " % (zygname))
        exit()

    try:
        fp3 = open(zygsort, "w")

    except (OSError, IOError) as e:
        print("ERROR: File Access Error!!!")
        print("ERROR:",e)
        print("Could NOT create File %s: " % (zygsort))
        exit()

    Files = os.listdir(os.getcwd())

    fileCnt = 0
    byteCnt = 0
    donor = {}
    for pname in Files:

        if pname.startswith("cluster-") and pname.endswith("-shared.vcf") and os.path.isfile(pname):
            pname = pname.strip()
            clustL= pname.split("_")
            clustN= clustL[0]      #Get Cluster ID Name
            clustS = volvox_scan_getSampleName(pname)     #Sample ID Name

            line  = pname
            line  = line.replace("cluster-","")
            line  = line.strip()
            word  = line.split("_")
            csize = int(word[1])

            if csize > MAXCLUSTPOP:
                #print("Exceeds Relevant Cluster Size Limits: ",pname)
                continue

            sample = volvox_scan_getSampleID(pname)
            print("Processing Sample: ", sample)
            fileCnt += 1
            byteCnt += os.path.getsize(pname)
            print("Opening unzipped File: ", pname, os.path.getsize(pname),"bytes")
            with open(pname,"r") as fp:

                gcat   = gcta = gccg = 0
                Gcat   = Gcta = Gccg = 0
                atgc   = atcg = atta = 0
                Atgc   = Atcg = Atta = 0
                others = 0
                Others = 0
                total  = 0
                hetcnt = 0
                homcnt = 0
                while True:
                    line = fp.readline()

                    if line == "":
                        break
                    
                    if line[0] == "#":
                        continue

                    line = line.strip()
                    word = line.split()

                    ref  = word[3].strip()
                    ale  = word[4].strip()
                    qual = float(word[5].strip())

                    if (qual < QUAL) or  (len(ref) != 1) or (len(ale) != 1) :
                        #print("Below Mutation Spectrum Criteria: ", qual, ref, ale)
                        continue
                    else:

                        if (ref == "G" and ale == "A") or (ref == "C" and ale == "T"):
                            gcat += 1

                        elif (ref == "G" and ale == "T") or (ref == "C" and ale == "A"):
                            gcta += 1

                        elif (ref == "G" and ale == "C") or (ref == "C" and ale == "G"):
                            gccg += 1

                        elif (ref == "A" and ale == "G") or (ref == "T" and ale == "C"):
                            atgc += 1

                        elif (ref == "A" and ale == "C") or (ref == "T" and ale == "G"):
                            atcg += 1

                        elif (ref == "A" and ale == "T") or (ref == "T" and ale == "A"):
                            atta += 1

                        else:
                            others += 1

                        total += 1
                        
                        if ("0/1" in line) or ("0|1" in line):
                            hetcnt += 1
                        elif ("1/1" in line) or ("1|1" in line):
                            homcnt += 1
                        else:
                            pass

            fp.close()

            if total > 0:
                Gcat   = gcat  *100/total
                Gcta   = gcta  *100/total
                Gccg   = gccg  *100/total
                Atgc   = atgc  *100/total
                Atcg   = atcg  *100/total
                Atta   = atta  *100/total
                Others = others*100/total

            fp1.write("%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" % \
                    (sample,Gcat,Gcta,Gccg,Atgc,Atcg,Atta,Others,\
                     gcat,gcta,gccg,atgc,atcg,atta,others,total))

            hcnt = hetcnt + homcnt

            if hcnt > 0:
                percent_het = 100*hetcnt/hcnt
            else:
                percent_het = -1
            fp2.write("%s\t%.2f\n" % (sample,percent_het))

            stuple = (clustS,percent_het)

            ans = donor.get(clustN,"None")

            if ans =="None":
                donor[clustN] = [stuple]
            else:
                ans.append(stuple)
                donor[clustN] = ans

        else:
            #print("Sneaky File or SubFolder: ", pname)
            pass

    print("\n\n")

    fp1.close()
    fp2.close()


    for k in donor:
        donor[k] = sorted(donor[k],key=lambda x:(x[1],x[0]))


    for key in donor:
        fp3.write("%s\t" %(key))
        for line in  donor[key]:
            x,y = line
            fp3.write("%s : %.1f%%\t" % (x,y))
        fp3.write("\n")

    fp3.close()

    return      #End of volvox_scan_specrum_zygosity

def volvox_scan_getParam(VALUE, promptStr):
    val = -1
    while val < 1:
        val_string = input(promptStr)

        if len(val_string) == 0:
            val = VALUE
            break
        else:

            try:

                val = float(val_string)
            except:
                print("ERROR: Invalid Input. Enter Positive Number: ", val_string)
                continue

    return val      #End of volvox_scan_getParam

def volvox_scan_getPromptStrings():
    pList = []
    print("\nFor Default Input Parameters ... Press the <Enter> Key\n")
    pList.append("Minimum Number of Shared Variants [250]: ")
    pList.append("Minimum Variant Quality-Score [10]: ")
    pList.append("Maximum No. of Samples Per Cluster  [5]: ")
    pList.append("Use SNPS Only <1> | All Variants <0> [1]: ")
    pList.append("Filter Common Variants <1> | Skip Filtering <0> [0]: ")
    
    return pList    #End of volvox_scan_getPromptStrings

def volvox_scan_getValues(VALUE,promptStr):
    val = -1
    while val < 0 or val > 1:
       
        val_string = input(promptStr)

        if len(val_string) == 0:
            val = VALUE
            break
        else:

            try:
                val = int(val_string)
            except:
                print("ERROR: Invalid Integer Value. Enter Zero or One.")
                continue

    return val      #End of volvox_scan_getValues
    
def volvox_scan_printOverFlow(overFlow,prefix):

    if len(overFlow) < 1:
        return

    pname = "Unused_Multi_Variants.txt"

    if len(prefix) > 0:
        pname = prefix.strip()+"-" + pname

    try:
        fp = open(pname,"w")

    except (OSError, IOError) as e:
        print("ERROR: File Access Error!!!")
        print("ERROR:",e)
        print("Could NOT create File %s: " % (pname))
        exit()

    for line in overFlow:
        line += "\n"
        fp.write(line)

    fp.close()

    return          #End of volvox_scan_printOverFlow

def volvox_scan_printStat(dataStat,mutStat,uniqCnt,overFlow,qual,prefix,Date):

    if len(dataStat) != 5:
        print("\n\nERROR: Inadmissible VCF File Analysis Results \n\n")
        exit()

    if len(mutStat) != 8:
        print("\n\nERROR: Invalid Cluster Analysis Results \n\n")
        exit()

    pname = "Analysis_Summary.txt"

    if len(prefix) > 0:
        pname = prefix.strip()+"-" + pname

    try:
        fp = open(pname,"w")

    except (OSError, IOError) as e:
        print("ERROR: File Access Error!!!")
        print("ERROR:",e)
        print("Could NOT create File %s: " % (pname))
        exit()

    print("########################################")
    print("#                                      #")
    print("#      Volvox Scan Summary Report      #")
    print("#                                      #")
    print("########################################")
    print("\nStart Time:",Date)

    #print("\n\nVolvox Scan Summary Report:")
    print("\n[File Processing]\n")
    print("Number of VCF Files: %d" % (dataStat[3]))
    print("Total File Size: %d Bytes" % (dataStat[4]))
    print("Total Number of Variants: %d" % (dataStat[0]))
    if dataStat[3] > 0:
        print("Average Number of Variants: %.2f" % (dataStat[0]/dataStat[3]))
    else:
        print("\nAverage Number of Variants: %.2f" % (0))
    #print("\n\nCluster Analysis:")
    print("Skipped MVSS Variants: %d" % (len(overFlow)))
    print("Total Number of SNPs: %d" %  (dataStat[1]))
    print("Total Number of INDELs: %d" % (dataStat[2]))
    print("Total Number of MNPs, OTHERs: %d" % (dataStat[0] -(dataStat[1] + dataStat[2])))

    fp.write("########################################\n")
    fp.write("#                                      #\n")
    fp.write("#      Volvox Scan Summary Report      #\n")
    fp.write("#                                      #\n")
    fp.write("########################################\n")
    fp.write("\nStart Time: %s\n" %(Date))
    fp.write("\n[File Processing]\n\n")
    fp.write("Number of VCF Files: %d\n" % (dataStat[3]))
    fp.write("Total File Sizes: %d Bytes\n" % (dataStat[4]))
    fp.write("Total Number of Variants: %d\n" % (dataStat[0]))
    if dataStat[3] > 0:
        fp.write("Average Number of Variants: %.2f\n" % (dataStat[0]/dataStat[3]))
    else:
        fp.write("Average Number of Variants: %.2f\n" % (0))
    fp.write("Skipped MVSS Variants: %d\n" % (len(overFlow)))
    fp.write("Total Number of SNPs: %d\n" %  (dataStat[1]))
    fp.write("Total Number of INDELs: %d\n" % (dataStat[2]))
    fp.write("Total Number of MNPs, OTHERs: %d\n" % (dataStat[0] -(dataStat[1] + dataStat[2])))
    
    if uniqCnt >= 0:
        print("Total Number of %s Unique Mutations: %d" % ("Q"+str(int(qual)),uniqCnt))
        fp.write("Total Number of %s Unique Mutations: %d\n" % ("Q"+str(int(qual)),uniqCnt))

        if dataStat[3] > 0:
            print("Average Number of %s Unique Mutations: %.2f" % ("Q"+str(int(qual)),uniqCnt/dataStat[3]))
            fp.write("Average Number of %s Unique Mutations: %.2f\n" % ("Q"+str(int(qual)),uniqCnt/dataStat[3]))
        else:
            print("Average Number of %s Unique Mutations: %.2f" % ("Q"+str(int(qual)),0))
            fp.write("Average Number of %s Unique Mutations: %.2f\n" % ("Q"+str(int(qual)),0))

    print("\n")
    print(" ---------------------")
    print("| Clustering Analysis |")
    print(" ---------------------")
    print("\n[All Clusters]\n")
    print("Total Number of Clusters: %d" % (mutStat[0]))
    print("Total Number of Shared Positions: %d" % (mutStat[2]))

    fp.write("\n\n")
    fp.write(" ---------------------\n")
    fp.write("| Clustering Analysis |\n")
    fp.write(" ---------------------\n")
    #fp.write("\n\nCluster Analysis]\n")
    fp.write("\n[All Clusters]\n\n")
    fp.write("Total Number of Clusters: %d\n" % (mutStat[0]))
    fp.write("Total Number of Shared Positions: %d\n" % (mutStat[2]))

    if mutStat[0] > 0:
        print("Average Number of Mutants in Cluster: %.2f" % (mutStat[1]/mutStat[0]))
        print("Average Number of Shared Mutations in Cluster: %.2f" % (mutStat[2]/mutStat[0]))

        fp.write("Average Number of Mutants in Cluster: %.2f\n" % (mutStat[1]/mutStat[0]))
        fp.write("Average Number of Shared Mutations in Cluster: %.2f\n" % (mutStat[2]/mutStat[0]))
    else:
        print("Average Number of Mutants in Cluster: %.2f" % (0))
        print("Average Number of Shared Mutations in Cluster: %.2f" % (0))

        fp.write("Average Number of Mutants in Cluster: %.2f\n" % (0))
        fp.write("Average Number of Shared Mutations in Cluster: %.2f\n" % (0))

    print("Total (Redundant) Number of Shared Mutations: %d" % (mutStat[3]))
    fp.write("Total (Redundant) Number of Shared Mutations: %d\n" % (mutStat[3]))

    print("\n[Retained Clusters]")
    print("\nTotal Number of Clusters: %d" % (mutStat[4]))
    print("Total Number of Shared Positions: %d" % (mutStat[6]))

    fp.write("\n[Retained Clusters]\n")
    fp.write("\nTotal Number of Clusters: %d\n" % (mutStat[4]))
    fp.write("Total Number of Shared Positions: %d\n" % (mutStat[6]))

    if mutStat[4] > 0:
        print("Average Number of Mutants in Cluster: %.2f" % (mutStat[5]/mutStat[4]))
        print("Average Number of Shared Mutations in Cluster: %.2f" % (mutStat[6]/mutStat[4]))

        fp.write("Average Number of Mutants in Cluster: %.2f\n" % (mutStat[5]/mutStat[4]))
        fp.write("Average Number of Shared Mutations in Cluster: %.2f\n" % (mutStat[6]/mutStat[4]))
    else:
        print("Average Number of Mutants in Cluster: %.2f" % (0))
        print("Average Number of Shared Mutations in Cluster: %.2f" % (0))

        fp.write("Average Number of Mutants in Cluster: %.2f\n" % (0))
        fp.write("Average Number of Shared Mutations in Cluster: %.2f\n" % (0))

    print("Total (Redundant) Number of Shared Mutations: %d" % (mutStat[7]))
    fp.write("Total (Redundant) Number of Shared Mutations: %d\n" % (mutStat[7]))
   
    Date = volvox_scan_printDate() 

    print("\nCompletion Time:",Date)
    fp.write("\nCompletion Time: %s\n" %(Date))

    print("\n\n")

    fp.close()

    return      #End of volvox_scan_printStat
    
def volvox_scan_filterVCFs(FILES,SNPS,counts,MINQUAL):

    uniqCnt = 0

    try:
        fp  = open(FILES,"r")
        print("List of VCF Files: ",FILES)

    except (OSError, IOError) as e:
        print("ERROR: File Access Error!!!")
        print("ERROR:",e)
        print("File Name %s: " % (FILES))
        exit()

    while True:
        line = fp.readline()

        if line =="":
            break

        FILE = line.strip()

        print("Filtering VCF File: ", FILE)

        ucnt = volvox_scan_filterVariants(FILE,SNPS,counts,MINQUAL)
        uniqCnt += ucnt

    return uniqCnt      #End of volvox_scan_filterVCFs

def volvox_scan_filterVariants(FILE,SNPS,counts,MINQUAL):

    tname = os.path.basename(FILE)

    if tname.endswith(".vcf.gz"):
        gzipped = True
        name    = tname.split(".vcf.gz")
        fname   = name[0]+"-uniq.vcf.gz"

    elif tname.endswith(".vcf"):
        gzipped = False
        name    = tname.split(".vcf")
        fname   = name[0]+"-uniq.vcf"

    else:
        print("\n\n***ERROR: Unknown File VCF File Name Extension" +tname+"\n\n")
        exit()
            
    sample = name[0].strip()

    if gzipped:
        try:
            fp  = gzip.open(FILE,"rt")
            fp1 = gzip.open(fname,"wb")
            print("Found a Gzipped File: ",FILE)

        except (OSError, IOError) as e:
            print("ERROR: File Access Error!!!")
            print("ERROR:",e)
            print("Files %s and %s: " % (FILE,fname))
            exit()
    else:
        try:
            fp  = open(FILE,"r")
            fp1 = open(fname,"w")
            print("Found an Unzipped File: ",FILE)

        except (OSError, IOError) as e:
            print("ERROR: File Access Error!!!")
            print("ERROR:",e)
            print("Files %s and %s: " % (FILE,fname))
            exit()

    uniqCnt = 0
    #uniqFile
    if gzipped:

        while True:

            line = fp.readline()

            if line == "":
                break

            if line[0] =="#":
                fp1.write(line.encode())
                continue
    
            word   = line.split("\t")
            posn   = word[0].strip()+":"+word[1].strip()
            ref    = word[3].strip()
            allele = word[4].strip()
            qual   = float(word[5].strip())

            if qual < MINQUAL:
                continue

            if (SNPS == 1) and  (len(ref) != 1 or len(allele) != 1):
                continue    #skip

            #Check Counts dictionary:
            result = counts.get(posn, "None")

            if result == "None":                # Found An OverFlow 
                fp1.write(line.encode())
                uniqCnt += 1
            else:
                if result == 1:
                    fp1.write(line.encode())
                    uniqCnt += 1
                else:
                    pass

    else:

        while True:

            line = fp.readline()

            if line == "":
                break

            if line[0] =="#":
                fp1.write(line)
                continue
    
            word   = line.split("\t")
            posn   = word[0].strip()+":"+word[1].strip()
            ref    = word[3].strip()
            allele = word[4].strip()
            qual   = float(word[5].strip())

            if qual < MINQUAL:
                continue

            if (SNPS == 1) and  (len(ref) != 1 or len(allele) != 1):
                continue    #skip

            #Check Counts dictionary:
            result = counts.get(posn, "None")

            if result == "None":                # Found An OverFlow 
                fp1.write(line)
                uniqCnt += 1
            else:
                if result == 1:
                    fp1.write(line)
                    uniqCnt += 1
                else:
                    pass

    fp.close()
    fp1.close()

    return uniqCnt      #End of volvox_scan_filterVariants

def volvox_scan_sortVCF():
    DIR=os.getcwd()
    print(DIR)
    Files = os.listdir(DIR)

    for File in Files:

        if File.endswith("-shared.vcf"):
            print("Sorting Shared VCF File Content: ", File)
            lines = []
            os.rename(File,"shared.txt")
            fp  = open("shared.txt","r")
            fp1 = open(File,"w")

            while True:
                line = fp.readline()
            
                if line == "":
                    break
    
                if line[0] == "#":
                    fp1.write(line)

                else:
                    line = line.split()
                    lines.append(line)

            slines= sorted(lines, key=lambda x: (x[0], int(x[1])))

            for line in slines:
                data = '\t'.join(line)
                fp1.write(data+"\n")

            fp.close()
            fp1.close()

    if os.path.exists("shared.txt"):
        os.remove("shared.txt")

    return      #End of volvox_scan_sortVCF

def volvox_scan_printDate():
    Day  = datetime.datetime.today().strftime("%A")
    Now  = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    Date = str(Day)+" " +str(Now)

    return Date     #End of volvox_scan_printDate

if __name__ == "__main__":
    main()




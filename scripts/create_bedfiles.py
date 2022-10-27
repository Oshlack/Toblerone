import sys
import itertools
import numpy as np
import csv

#todo argparse to clean this up   
inputbed12 = sys.argv[1]
outputdir = sys.argv[2]

includetranscripts=[]
if len(sys.argv) >= 4:
    includefile = sys.argv[3]
    subset=True
    with open(includefile, 'r') as include_file_handle:
        includetranscripts = include_file_handle.readlines()

else:
    subset=False



#print(subset)
#quit()




full_gene_dict = {}
with open(inputbed12, 'r') as csv_file:
    for row in csv.reader(csv_file, delimiter='\t'):
        full_gene_dict[row[3]] = '\t'.join(row)
        #break
        
def get_internals(combo):
    #create an index of the indices to remove, e.g the dels
    subs=[]
    for i in range(1,len(combo)):   # for the length of the full sequence
        for y in range(1,i):
            #print(list(range(y,i)))
            # starting at 1, the second, up to i
            mask  =  np.array(combo[range(y,i)])  # get indices of that range
            #print(mask)
            subs.append(combo[np.in1d(combo, mask, invert=True)])
    
    return(subs)



def make_bed_file(lengths, positions, orig_exon, seed_bed1, seed_bed2, bedfilename, gapsallowed = False):
    
    if (not gapsallowed):
        bedfilename_mod = bedfilename+"_nogaps_optimized"
    else:
        bedfilename_mod = bedfilename+"_all_new"
   
    lengths_lookup = dict(zip(range(1,len(lengths)+1),lengths))
    positions_lookup = dict(zip(range(1,len(lengths)+1),positions))

 
    all_subsets_v2 = get_internals(np.arange(1,1+len(lengths_lookup)))

    total_written = 0
    #print(all_subsets)



    for sub in all_subsets_v2:
        hits = sub
        name = "_".join(  [str(x) for x in hits])
        del_name = [str(x) for x in range(1,1+len(lengths_lookup)) if x not in hits]

        a_actual = [lengths_lookup[int(i)] for i in sub]
        b_actual = [positions_lookup[int(i)] for i in sub]

        del_name = "_".join(   del_name)

        if(not name):
            name = "error_no_name"
            sys.stderr.write("No del found: ",seed_bed1,seed_bed2)
        if(del_name): # no name means no del, so all internal exons skipped
            #print(seed_bed1+"del"+del_name+"\t"+seed_bed2)
            print(seed_bed1+"_del"+del_name+"\t"+seed_bed2+"\t"+str(len(list(a_actual)))+"\t"+",".join([str(s) for s in a_actual])+","+"\t"+",".join([str(s) for s in b_actual])+",",file=open(outputdir+"/"+bedfilename_mod+".bed", "a"))
            total_written = total_written +1
        else:
            print(seed_bed1+"_del"+del_name+"\t"+seed_bed2+"\t"+str(len(list(a_actual)))+"\t"+",".join([str(s) for s in a_actual])+","+"\t"+",".join([str(s) for s in b_actual])+",",file=open(outputdir+"/"+bedfilename_mod+".skipped", "a"))
            sys.stderr.write("No del found: ",seed_bed1,seed_bed2)

#print(includetranscripts)
#print(len(full_gene_dict.items()))
for k,v  in full_gene_dict.items():
    bedfilename = k
    #print(k)
    #print(v)


    # Add canonical transcrit bed without the _del suffix
    print(v,file=open(outputdir+"/"+k+"_nogaps_optimized"+".bed", "a"))
    bed_contents = v.split("\t")


    # If Subset mode, only entreis in the bed file that match will have there values generated
    # transcripts to match column 4
    if subset:
     #   print(bed_contents[3])
      #  print(includetranscripts[0])
       # print(any(s in bed_contents[3]  for s in includetranscripts))
        if not any(s.strip("\n") in bed_contents[3]  for s in includetranscripts):
            #print("NOT IN")
            continue
    #print(bed_contents)


    #continue
    seed_bed1 = "\t".join(bed_contents[0:4])
    seed_bed2 = "\t".join(bed_contents[4:9])
    orig_exon = bed_contents[9]
    lengths_s = bed_contents[10]
    positions_s = bed_contents[11]
    lengths = [int(s) for s in lengths_s.split(',')[0:int(orig_exon)]]
    positions = [int(s) for s in  positions_s.split(',')[0:int(orig_exon)]]

    make_bed_file(lengths, positions, orig_exon, seed_bed1, seed_bed2, bedfilename, gapsallowed = False)
    #make_bed_file(lengths, positions, orig_exon, seed_bed1, seed_bed2, bedfilename, gapsallowed = True)



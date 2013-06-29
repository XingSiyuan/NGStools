import argparse
import sys
import os
import re
import mysql.connector

parser = argparse.ArgumentParser( description='Aligns fasta file and returns stats per aligned base')
parser.add_argument("-p", "--positionfile", help="file with 60K SNP positions on build 10.2", nargs=1)
parser.add_argument("-o", "--outfile", help="output file", nargs=1)

db = mysql.connector.Connect(user='xxxxx',password='xxxxx',host='localhost',database='pig_hapmap2')
cursor = db.cursor()

def select_genotypes(snp_name):
   output={}
   #print(snp_name)
   stmt_select = "select dna_name,fwallele1,fwallele2 from allgenotypes16 where SNP = '"+snp_name+"'"
   cursor.execute(stmt_select)

   for row in cursor.fetchall():
      output[row[0]]=row[1]+row[2]

   return output


def select_ids():
   output={}
   stmt_select = "SELECT dna_name from sample_sheet8 where callrate > 0.7"
   cursor.execute(stmt_select)

   for row in cursor.fetchall():
      output[row[0]]=0

   return output

def read_positions (filename):
   pos={}
   with open(filename) as file:
      for  l in file.readlines():
         content = l.split()
         #print(content)
         pos[content[1]+'_'+content[2]]=content[0]

   return pos

args = parser.parse_args()
pos_file=args.positionfile[0]
alignedfile=args.outfile[0]

kpos = read_positions(pos_file)

iddict = select_ids()
print(len(kpos))
print(len(iddict))
totalsnps=0
for l in sys.stdin.readlines():
   l=l.rstrip().split()
   chrom=l[0]
   refpos=l[1]
   refbase=l[2]
   consensusbase=l[3]
   consensusqual=l[4]
   snpqual=l[5]
   mappingqal=l[6]
   readdepth=l[7]
   readbases=l[8]
   basesqual=l[9]
   if int(snpqual) > 20 and int(consensusqual) > 20 and re.match('[RYSWKM]',consensusbase):
      if chrom+'_'+refpos in kpos:
         snp=kpos[chrom+'_'+refpos] 
         genotypes=select_genotypes(snp)
         if len(genotypes)>0:
            totalsnps+=1
            for ind in genotypes.keys():
               if ind in iddict:
                  alleles = list(genotypes[ind])
                  if alleles[0] != alleles[1]:
                     iddict[ind]+=1 

for ind in iddict.keys():
   pct = int(iddict[ind])/totalsnps
   if pct > 0.7:
       print(ind,pct,sep='\t')

cursor.close()
db.close()

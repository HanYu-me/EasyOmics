import sys
file1=sys.argv[1]
file2=sys.argv[2]
f_in=open(file1)
f_out=open(file2,"w")
a=f_in.readline()
while a:
    a=a.split("\t")
    a[0]=a[0].replace("Chr","")
    a[0]=a[0].replace("chr","")
    a[0]=a[0].replace("Chromosome","")
    a[0]=a[0].replace("chromosome","")
    a[0]=a[0].replace("Ay","")
    a[0]=a[0].replace("_RagTag","")
    a='\t'.join(a)          
    f_out.writelines(a)
    a=f_in.readline()
f_in.close()
f_out.close()

import sys
file1=sys.argv[1]
file2=sys.argv[2]
f_in=open(file1)
f_out=open(file2,"w")
a=f_in.readline()
while a:
    
    if a[0]!= "#":
        if "scaffold" in a:
            a=f_in.readline()
            continue
        if "*" in a:
            a=f_in.readline()
            continue
        a=a.split("\t")
        sub=a[2].split(":")
        a[1]=sub[1]
        a='\t'.join(a)  
    f_out.writelines(a)
    a=f_in.readline()
f_in.close()
f_out.close()

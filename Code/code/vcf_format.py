import sys
file1=sys.argv[1]
file2=sys.argv[2]
f_in=open(file1)
f_out=open(file2,"w")
a=f_in.readline()
while a:
    if a[0]=="#" and a[1]!="#":
        a=a.replace("|","_")
        a=a.split("\t")
        for i in range(len(a)):
            if i>=9:
                if "_" not in a[i]:
                    a[i]=a[i]+"_"+a[i]
        a='\t'.join(a)
        #a=a+"\n"
    
    if a[0]!= "#":
        if "scaffold" in a:
            a=f_in.readline()
            continue
        if "*" in a:
            a=f_in.readline()
            continue
        
        a=a.replace("Chr","")
        a=a.replace("chr","")
        a=a.replace("Chromosome","")
        a=a.replace("chromosome","")
        a=a.split("\t")
        if "rs" not in a[2] and ":" not in a[2]:
            a[2]=a[0]+":"+a[1]
        a='\t'.join(a)  
        #a=a+"\n"
        
    f_out.writelines(a)
    a=f_in.readline()
f_in.close()
f_out.close()

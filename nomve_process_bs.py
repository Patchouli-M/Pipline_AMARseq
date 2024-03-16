from multiprocessing import Pool
import pandas as pd 
import sys
import os
import time

# need dirs in sample_lib/data named sample_name
# and sample_1.fq.gz in sample_2.fq.gz in sample_name

tool_lib=sys.path[0]
sample_lib="~/MyFiles/20240203"
reference_lib="~/MyFiles/reference/bismark_reference"
USE_CORE = 35
def get_sample(dir:str=sample_lib):
    data_dir=os.path.join( dir, "data" )
    print ("=========== Check Sample ===========")
    global sample_l
    sample_l=[]
    for root,dir,file in os.walk(data_dir):
        if (not file):
            continue
        file.sort()
        if "MD5.txt" in file :
            file.remove('MD5.txt')
        if (len(file)!=2):
            print ("Error sample, check please. ")
            sys.exit()
        print(file)
        for i in range(len(file)):
            if ( (file[i].split('_')[-1]) != str(i+1)+".fq.gz" ):
                print ("Error sample, check please. ")
                sys.exit()
        f1 = file[0].split('.fq.gz')[0][:-2]
        f2 = file[1].split('.fq.gz')[0][:-2]
        if (f1 != f2 ):
            print ("Error sample, check please. ")
            sys.exit()
        sample_l.append(f1)
        print (file)
    sample_l.sort()
    print(sample_l)

def multi_sys_cmd(cmd:str):
    print(cmd)
    os.system(cmd)

def trim():
    # Input: Raw data
    # Output : trimmed data 
    print ("=========== Trim ===========")
    out_dir=os.path.join(sample_lib,"trimmed")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    cmd_l=[]
    for f in sample_l:
        f1 = f+"_1.fq.gz"
        f2 = f+"_2.fq.gz"
        f1 = os.path.join(sample_lib,"data",f,f1)
        f2 = os.path.join(sample_lib,"data",f,f2)
        cmd="trim_galore --quality 20 --phred33 --stringency 3 --gzip --length 36 --paired --cores 4 --output_dir %s %s %s " % (out_dir,f1,f2)
        cmd_l.append(cmd)
    core_num = USE_CORE
    pool = Pool(core_num)
    pool.map(multi_sys_cmd, cmd_l)
    pool.close()
    print ("======== Trim Finish ========")

def bismark():
    # Input: trimmed data 
    # Output : bismark bam
    print ("=========== bismark compare ===========")
    out_dir=os.path.join(sample_lib,"bismark")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    cmd_l=[]
    for f in sample_l:
        f1 = f+"_1_val_1.fq.gz"
        f2 = f+"_2_val_2.fq.gz"
        f1 = os.path.join(sample_lib,"trimmed",f1)
        f2 = os.path.join(sample_lib,"trimmed",f2)
        cmd="bismark -fastq --output_dir  %s --temp_dir %s/tmp --multicore 8 --non_directional -bowtie2 %s -1 %s -2 %s" % (out_dir,out_dir,reference_lib,f1,f2)
        cmd_l.append(cmd)
    # pool = Pool(USE_CORE)
    # pool.map(multi_sys_cmd, cmd_l)
    # pool.close()
    for cmd_str in cmd_l:
        multi_sys_cmd(cmd_str)
    print ("======== bismark compare Finish ========")

def bismark_dup():
    work_dir=os.path.join(sample_lib,"bismark")
    raw_dir=os.path.join(sample_lib,"bismark_raw")
    if not os.path.exists(raw_dir):
        os.makedirs(raw_dir)
    cmd = 'cp ' + work_dir + '/* ' + raw_dir
    print (cmd)
    os.system(cmd)
    print ("======== cp Raw Finish ========")
    for f in sample_l:
        f1 = f+"_1_val_1.fq.gz_bismark_bt2_pe.bam"
        f1 = os.path.join(work_dir,f1)
        cmd = 'deduplicate_bismark -p --bam '+f1
        print(cmd)
        os.system(cmd)
        d_f = f+"_1_val_1.fq.gz_bismark_bt2_pe.deduplicated.bam"
        d_f = os.path.join(work_dir,d_f)
        cmd = 'mv '+ d_f + ' ' +f1
        print (cmd)
        os.system(cmd)

def extractor_C():
    # Input: bismark bam
    # Output : covs of CG and all_C
    print ("=========== extractor_C ===========")
    out_dir_allC=os.path.join(sample_lib,"extractor_all_C")
    if not os.path.exists(out_dir_allC):
        os.makedirs(out_dir_allC)
    out_dir_CpG=os.path.join(sample_lib,"extractor_CpG")
    if not os.path.exists(out_dir_CpG):
        os.makedirs(out_dir_CpG)
    
    cmd_CpG_l=[]
    cmd_all_C_l=[]
    cmd_c2c_l=[]
    cmd_unzip_l=[]
    for f in sample_l:
        f1 = f+"_1_val_1.fq.gz_bismark_bt2_pe.bam"
        f1 = os.path.join(sample_lib,"bismark",f1)

        cmd_all_C = "bismark_methylation_extractor -p --buffer_size 128G --ample_memory --multicore 2 --output %s --cytosine_report --CX --CX --genome_folder %s %s" % (out_dir_allC,reference_lib,f1)

        cmd_CpG = "bismark_methylation_extractor -p --buffer_size 128G --ample_memory --multicore 2 --output %s --cytosine_report --genome_folder %s %s" % (out_dir_CpG,reference_lib,f1)
        
        cov_f=f+"_1_val_1.fq.gz_bismark_bt2_pe.bismark.cov.gz"
        cov_f = os.path.join(out_dir_allC,cov_f)
        
        # out_dir_allC 
        cmd_c2c="coverage2cytosine --genome_folder %s --dir %s -o %s --gc %s" %(reference_lib,out_dir_allC,f,cov_f)
        
        CpG_cov_f=f+"_1_val_1.fq.gz_bismark_bt2_pe.bismark.cov.gz"
        CpG_cov_f = os.path.join(out_dir_CpG,CpG_cov_f)

        cmd_unzip="gunzip %s" % (CpG_cov_f)

        cmd_CpG_l.append(cmd_CpG)
        cmd_unzip_l.append(cmd_unzip)
        cmd_all_C_l.append(cmd_all_C)
        cmd_c2c_l.append(cmd_c2c)
    core_num = USE_CORE
    pool = Pool(core_num)
    pool.map(multi_sys_cmd, cmd_CpG_l)
    pool.map(multi_sys_cmd, cmd_unzip_l)
    pool.map(multi_sys_cmd, cmd_all_C_l)
    pool.map(multi_sys_cmd, cmd_c2c_l)
    pool.close()
    print ("======== extractor_C Finish ========")

def multi_remove_GCG( obj_cov_f:str ):
        out_dir=os.path.join(sample_lib,"dropped")
        omGCG_f="omGCG.bed"
        GCG=pd.read_csv(os.path.join(tool_lib,omGCG_f),sep='\t',header=None,dtype={0:str,1:int,2:int})
        obj_cov=pd.read_csv(obj_cov_f,sep='\t',header=None,dtype={0:str,1:int,2:int,4:int,5:int})
        new_cov=pd.DataFrame()
        count=0
        chr_list=obj_cov.drop_duplicates([0])[0].to_list()
        chr_list.sort()
        
        for i in chr_list:
            print (i)
            temp_GCG = GCG[GCG[0]==i]
            temp_obj = obj_cov[obj_cov[0]==i]
            GCG_l=(temp_GCG[2]-1).to_list()
            new_cov=new_cov.append(temp_obj[~temp_obj[1].isin(GCG_l)])
        new_cov[1]=new_cov[1].astype(int)  
        new_cov[2]=new_cov[2].astype(int)  
        new_cov[4]=new_cov[4].astype(int)  
        new_cov[5]=new_cov[5].astype(int)
        print (new_cov)
        f_c=obj_cov_f.split('/')[-1]
        if ("GpC" in f_c):
            f_c=f_c.split(".GpC.cov")[0]
            f_c_o=f_c+".GpC.cov"
            f_c=f_c+".GpC.dropped.cov"
        else:
            f_c=f_c.split("_1_val_1.fq.gz_bismark_bt2_pe.bismark.cov")[0]
            f_c_o=f_c+".CpG.cov"
            f_c=f_c+".CpG.dropped.cov"
        f_c=os.path.join(out_dir,f_c)
        f_c_o=os.path.join(out_dir,f_c_o)
        new_cov.to_csv(f_c,sep='\t',header=False,float_format='%.6f',index=False)
        obj_cov.to_csv(f_c_o,sep='\t',header=False,float_format='%.6f',index=False)

def remove_GCG():
    print ("=========== Drop GCG ===========")
    out_dir=os.path.join(sample_lib,"dropped")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    drop_list=[]    
    for i in sample_l:
        CG_cov = i+"_1_val_1.fq.gz_bismark_bt2_pe.bismark.cov"
        CG_cov = os.path.join(sample_lib,"extractor_CpG",CG_cov)
        GC_cov = i+".GpC.cov"
        GC_cov = os.path.join(sample_lib,"extractor_all_C",GC_cov)
        drop_list.append(CG_cov)
        drop_list.append(GC_cov)
    core_num = USE_CORE
    print (drop_list)
    pool = Pool(core_num)
    pool.map(multi_remove_GCG, drop_list)
    pool.close()
    print ("=========== Drop GCG Finish ===========")

def trans_chr():
    print ("=========== transfor chr ===========")
    out_dir=os.path.join(sample_lib,"transed")
    inp_dir=os.path.join(sample_lib,"dropped")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    for root,dir,file in os.walk(inp_dir):
        for i in file:
            in_f=os.path.join(inp_dir,i)
            out_f=os.path.join(out_dir,i)
            cmd="grep '^[0-9]' %s > %s" % (in_f,out_f)
            multi_sys_cmd(cmd)
            cmd = "grep 'X' %s >> %s" % (in_f,out_f)
            multi_sys_cmd(cmd)
            cmd = "grep 'Y' %s >> %s" % (in_f,out_f)
            multi_sys_cmd(cmd)
            cmd = "sed -i -e 's/^/chr/' %s" % (out_f)
            multi_sys_cmd(cmd)
            cmd = "wc -l %s" % (in_f)
            multi_sys_cmd(cmd)
            cmd = "wc -l %s" % (out_f)
            multi_sys_cmd(cmd)
    print ("=========== transfor chr finish ===========")
    


get_sample() 
trim()
bismark() 
bismark_dup()
extractor_C()
remove_GCG() 
trans_chr()
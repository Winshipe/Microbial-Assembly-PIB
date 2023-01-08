configfile: "config.yaml"

import pandas
import subprocess
import re
import platform
import datetime

now = datetime.datetime.now()
time = now.strftime("%d-%m.%H%M")

"""Here is a pipeline to run all the steps needed to get to an iCAMP result. The starting point is a file which contains the following columns: ASV, ASV_SEQUENCE, DOMAIN, Phylum, Class, Order, Family, Genus, <SITE>_<DEPTH>_<NA type>..."""

#print("Starting iCAMP pipeline")

rule all:
    input:
        expand([\
            "{sample_id}.phylo.txt",\
            "{sample_id}.counts.txt",\
            "{sample_id}.nwk",\
            "{sample_id}.treat2col.txt",\
            "{sample_id}.env.txt",\
            "{sample_id}.selection.png",\
            "{sample_id}.tax.txt"\
	 ],sample_id = config["sample id"] + "-"+time)

rule filter_tax:
    input:
        "Lilian_ASV_table_taxonomy.txt"
    output:
        "{sample_id}.tax.txt"
    run:
        #print("FT")
        raw_tax = pandas.read_csv(input[0], sep="\t")
        cols = list(raw_tax.columns)
        sublist = [c for c in cols if "DNA" in c and "M5_A" in c and "DNA" in c]
        tax = raw_tax.loc[(raw_tax[sublist].sum(axis=1) >= 30) & (raw_tax["DOMAIN"] == "Archaea")]
        tax.to_csv(output[0],sep="\t",index=False)

rule make_tree:
# generate the required tree/newick file with the 16S information from our ASVs
    input:
        "{sample_id}.tax.txt"
    output:
        "{sample_id}.nwk"
    run:
        #print("MT")
        fas_path = "temp.fa"
        aln_path = "temp.aln.fa"
        tax = pandas.read_csv(input[0],sep="\t")

        with open(fas_path,"w") as f:
            for idx,r in tax.iterrows():
                f.write(">{ASV}\n{ASV_SEQUENCE}\n".format(**r))
        pfm = platform.system()
        #print(pfm)
        ft = "fasttree" if pfm != "Windows" else "FastTree.exe"
        msl = "muscle" if pfm != "Windows" else "muscle5.1.win64.exe"
        subprocess.run([msl,"-super5",fas_path, "-output",aln_path],stdout=subprocess.PIPE)
        subprocess.run([ft,"-out",output[0],"-nt","-gtr",aln_path],stdout=subprocess.PIPE)

rule comm:
# makes the comm file // gets the taxonomy information for the relevant ASVs
    input:
        "{sample_id}.tax.txt"
    output:
        "{sample_id}.phylo.txt"
    run:
        #print("C")
        tax = pandas.read_csv(input[0],sep="\t")
        cols = ["ASV","DOMAIN","Phylum","Class","Order","Family","Genus"]
        tax[cols].to_csv(output[0],sep="\t",index=False)
rule otus:
# make the otus file // gets ASV and counts
    input:
        "{sample_id}.tax.txt"
    output:
        "{sample_id}.counts.txt"
    run:
        #print("O")
        tax = pandas.read_csv(input[0],sep="\t")
        col_names = list(tax.columns)
        #site = config["site"]
        na_type = config["nucleic acid"]
        my_cols = [c for c in col_names if na_type in c]# (site in c) and (na_type in c)] # eg M5_A and DNA -> M5_A_5_DNA
        my_cols.insert(0,"ASV")
        tax[my_cols].to_csv(output[0],sep="\t",index=False)
rule make_env:
    input:
        "{sample_id}.tax.txt"
    output:
        environment = "{sample_id}.env.txt",
        treat2col = "{sample_id}.treat2col.txt"
    run:
        #print("ME")
        tax = pandas.read_csv(input[0],sep="\t")
        #site = config["site"]
        na_type = config["nucleic acid"]
        cutoff = float(config["depth cutoff"])
        depths = []
        rgx = re.compile("(?<=_)[0-9]+\.?[0-9]*(?=_)")
        my_cols = []
        col_names = list(tax.columns)
        for c in col_names:
            if na_type not in c: # eg M5_A and DNA -> M5_A_5x5_DNA
                continue
            c_ = c.replace("x",".")
            m = rgx.search(c_)
            if m:
                my_cols.append(c)
                depths.append(m.group(0))
        with open(output.environment,"w") as f:
            f.write("SampleID\tDepths\n")
            for s_id, depth in zip(my_cols, depths):
                f.write("{}\t{}\n".format(s_id, depth))

        with open(output.treat2col, "w") as f:
            f.write("SampleID\tManagement\tLocation\n")
            for s_id, depth in zip(my_cols, depths):
                site = s_id.split("_")[0]
                mgmt = "BT" if float(depth) > 6 else "NB"
                loc = site + "_" + mgmt
                f.write("{}\t{}\t{}\n".format(s_id, mgmt, loc))

rule iCAMP:
    input:
        clas = "{sample_id}.phylo.txt",
        comm = "{sample_id}.counts.txt",
        tree = "{sample_id}.nwk",
        t2c = "{sample_id}.treat2col.txt",
        envm = "{sample_id}.env.txt"
    params:
        pfx = "{sample_id}",
        ds = 0.2,
        binsize = 12
    output:
        "{sample_id}.selection.png"
    script:
        "run_iCAMP.R"

import os, subprocess
import numpy as np, pandas as pd
from pyfaidx import Fasta
import click
from pathos import multiprocessing as mp
import glob
from itertools import combinations
import re
import pickle

#################
#### CLASSES ####

#########################
#### SYNTENY CLASSES ####
class PairwiseSynteny:
    def __init__(self, q_genome, s_genome, synteny_file='', loci_threshold = 4):
        self.synteny_file = synteny_file
        self.q_genome = q_genome
        self.s_genome = s_genome
        self.loci_threshold = loci_threshold

    def generate_synteny_structure(self,synteny_path):
        """Take anchor file or synteny file and searches for starting and ending genes for each syntenic block"""
        if self.synteny_file.endswith('.unout'):
            self.unout2structure(self.q_genome, self.s_genome)
        elif self.synteny_file.endswith('.lifted.anchors'):
            self.anchor2structure(self.q_genome, self.s_genome)
        elif self.synteny_file.endswith('.bed'):
            self.import_synteny_structure()
        else:
            self.run_synteny(self.q_genome,self.s_genome,synteny_path)
            self.anchor2structure(self.q_genome, self.s_genome)
        self.synteny_structure_index()

    def synteny_structure_index(self):
        self.synteny_structure = self.synteny_structure.rename(dict(zip(range(self.synteny_structure.shape[0]),self.synteny_structure['q_chr']+'\t'+self.synteny_structure['q_xi']+'\t'+self.synteny_structure['q_xf'])))

    def import_synteny_structure(self):
        self.synteny_structure = pd.read_table(self.synteny_file,header=None,names=['q_chr','q_xi','q_xf','s_chr','s_xi','s_xf'])
        self.synteny_structure_index()
        #self.synteny_structure = self.synteny_structure.rename(dict(enumerate((self.synteny_structure['q_chr']+'\t'+self.synteny_structure['q_xi']+'\t'+self.synteny_structure['q_xf']).as_matrix().tolist())))

    def unout2structure(self,q_genome,s_genome):
        with open(self.synteny_file,'r') as f:
            lines = np.array(f.read().splitlines())
        anchors = np.array_split(lines,np.where(np.vectorize(lambda line: line.startswith('\t') == 0)(lines))[0])
        synteny_structure = []
        for anchor in anchors:
            anchor = np.array(map(lambda line: line.split('\t')[1].split(','),anchor[1:].tolist()))
            if len(anchor) >= self.loci_threshold and anchor.tolist():
                q_genes, s_genes = anchor[:,2], anchor[:,5]
                q_coords, s_coords = q_genome.df.loc[q_genes,:], s_genome.df.loc[s_genes,:]
                synteny_structure.append([q_coords.iloc[0,0],q_coords[['xi','xf']].values.min(),q_coords[['xi','xf']].values.max(),s_coords.iloc[0,0],s_coords[['xi','xf']].values.min(),s_coords[['xi','xf']].values.max()])
        self.synteny_structure = pd.DataFrame(synteny_structure,columns=['q_chr','q_xi','q_xf','s_chr','s_xi','s_xf']).astype(str)#,index = np.vectorize(lambda x: '\t'.join(map(str,x[:3])))(synteny_structure))

    def run_synteny(self,genome1,genome2, synteny_path):
        pwd = os.getcwd()
        os.chdir(synteny_path)
        #subprocess.call('rm {0}/*.bck {0}/*.prj {0}/*.sds {0}/*.ssp {0}/*.suf {0}/*.tis {0}/*.des {0}/*.bed {0}/*.cds'.format(synteny_path),shell=True)
        for abs_path, link_name in zip([genome1.bed_file,genome2.bed_file,genome1.CDS_file,genome2.CDS_file],[genome1.protID+'.bed',genome2.protID+'.bed',genome1.protID+'.cds',genome2.protID+'.cds']):
            subprocess.call('ln -s %s %s'%(abs_path,link_name),shell=True)
        try:
            subprocess.call('python -m jcvi.compara.catalog ortholog --no_strip_names %s %s'%(genome1.short_name,genome2.short_name),shell=True)
        except:
            subprocess.call('python -m jcvi.compara.catalog ortholog --no_strip_names %s %s'%(genome1.short_name,genome2.short_name),shell=True)
        if genome1.short_name != genome1.protID and genome2.short_name != genome2.protID:
            subprocess.call('mv %s.%s.lifted.anchors %s.%s.lifted.anchors'%(genome1.short_name,genome2.short_name,genome1.protID,genome2.protID),shell=True)
        self.synteny_file = os.path.abspath('%s.%s.lifted.anchors'%(genome1.protID,genome2.protID))
        os.chdir(pwd)

    def anchor2structure(self, q_genome, s_genome):
        anchor_file = self.synteny_file
        with open(anchor_file,'r') as f:
            anchors = f.read().split('###')[1:]
        synteny_structure = []
        for anchor in anchors:
            if anchor:
                genes = np.array([line.split()[:2] for line in anchor.splitlines() if line])
                if genes.shape[0] >= self.loci_threshold:
                    q_genes, s_genes = genes[:,0] , genes[:,1]
                    q_coords, s_coords = q_genome.df.loc[q_genes,:], s_genome.df.loc[s_genes,:]
                    #print q_coords[['xi','xf']]
                    synteny_structure.append([q_coords.iloc[0,0],q_coords[['xi','xf']].values.min(),q_coords[['xi','xf']].values.max(),s_coords.iloc[0,0],s_coords[['xi','xf']].values.min(),s_coords[['xi','xf']].values.max()])
        self.synteny_structure = pd.DataFrame(synteny_structure,columns=['q_chr','q_xi','q_xf','s_chr','s_xi','s_xf']).astype(str)

    def synteny_structure_2_bed(self,filename):
        self.synteny_structure.to_csv(filename,sep='\t',index=False,header=None)

    def synteny_structure_2_link(self, filename):
        df = self.synteny_structure
        df['q_chr'] = np.vectorize(lambda x: self.q_genome.protID+'-'+x)(df['q_chr'])
        df['s_chr'] = np.vectorize(lambda x: self.s_genome.protID+'-'+x)(df['s_chr'])
        df.to_csv(filename,sep=' ',index=False,header=None)
        self.link = os.path.abspath(filename)

    def export_karyotypes(self,circos_input):
        self.s_genome.export_karyotype(circos_input+'/'+self.s_genome.protID+'.karyotype.txt')
        self.q_genome.export_karyotype(circos_input+'/'+self.q_genome.protID+'.karyotype.txt')

######################
#### GENOME CLASS ####

class Genome:
    def __init__(self, fasta_file, bed_file, protID, gff_file = '', gene_info = 'Name'):
        self.fasta_file = fasta_file
        self.fasta = Fasta(fasta_file)
        self.gene_info = gene_info
        self.bed_file = os.path.abspath(bed_file)
        self.short_name = self.bed_file.split('/')[-1].replace('.bed3','').replace('.bed','')
        self.protID = protID
        self.CDS_file = self.bed_file.replace('.bed3','.cds').replace('.bed','.cds')
        self.gff_file = gff_file
        if self.gff_file and os.path.exists(self.bed_file) == 0 or (os.path.exists(self.bed_file) and os.stat(self.bed_file).st_size == 0):
            #click.echo('python -m jcvi.formats.gff bed --type=mRNA --key=%s %s > %s'%(self.gene_info,self.gff_file,self.bed_file))
            subprocess.call('python -m jcvi.formats.gff bed --type=mRNA --key=%s %s -o %s'%(self.gene_info,self.gff_file,self.bed_file),shell=True)
            """
            with open(gff_file,'r') as f:
                for header_line,line in enumerate(f):
                    if line.startswith('#') == 0:
                        break
            df = pd.read_table(gff_file, skiprows=header_line, header=None,names=['chr','rm_1','feature','xi','xf','rm_3','rm_4','rm_5','Annotation'])
            df = df[df['feature'] == 'mRNA'].drop([feat for feat in list(df) if 'rm' in feat],axis=1)
            df = df[np.vectorize(lambda line: 'longest=1' in line)(df['Annotation']).astype(bool)]
            df['Gene'] = np.vectorize(lambda line: line.split(';')[1].replace('Name=',''))(df['Annotation'])
            df['xi'] -= 1
            df = df.drop(['feature','Annotation'],axis=1).reindex(columns=['chr','xi','xf','Gene'])
            self.df = df
            """
        self.df = pd.read_table(self.bed_file,header=None,names=['chr','xi','xf','Gene'],dtype={'chr':str,'xi':np.int,'xf':np.int,'Gene':str},usecols=[0,1,2,3])
        self.df = self.df.set_index('Gene')

    def export_bed(self,filename):
        df = self.df.reset_index().rename(dict(index='Gene'),axis='columns').reindex(columns=['chr','xi','xf','Gene'])
        df.to_csv(filename,sep='\t',index=False,header=None)

    def extract_CDS(self): # python -m jcvi.formats.gff uniq t.PAC2_0.316.gff3 -o uniq.gff3
        if os.path.exists(self.CDS_file) == 0 or (os.path.exists(self.CDS_file) and os.stat(self.CDS_file).st_size == 0):
            subprocess.call('python -m jcvi.formats.gff load %s %s --parents=mRNA --children=CDS --id_attribute=%s -o %s'%(self.gff_file,self.fasta_file,self.gene_info,self.CDS_file),shell=True)

    def export_karyotype(self, filename, n_chromosomes=25, shorten_chr=False):
        df = pd.read_table(self.fasta_file+'.fai', header=None,names=['chr','length'],usecols=[0,1],dtype=dict(zip(['chr','length'],[str,np.int])))
        df = df.sort_values(['length'],ascending=False)
        if n_chromosomes < df.shape[0]:
            df = df.iloc[:n_chromosomes,:].reset_index(drop=True)
        out_txt = []
        for i in range(df.shape[0]):
            chrom = df.loc[i,'chr']
            chr_name = chrom if not shorten_chr else chrom[0] + chrom.split('_')[-1]
            if i >= 25:
                out_txt.append('chr - %s-%s %s 0 %d %d,%d,%d\n'%(self.protID,chrom,chr_name,df.loc[i,'length'],randint(1,255),randint(1,255),randint(1,255)))
            else:
                out_txt.append('chr - %s-%s %s 0 %d chr%d\n'%(self.protID,chrom,chr_name,df.loc[i,'length'],i+1))
        with open(filename,'w') as f:
            f.writelines(out_txt)
        self.karyotype = os.path.abspath(filename)

######################
#### CIRCOS CLASS ####

class Circos:
    def __init__(self,PairwiseSynteny):
        self.synteny = PairwiseSynteny

    def write_ideogram_config(self, filename='txideogram.conf'):
        with open(filename,'w') as f:
            f.write("""<ideogram>
                show = yes
                <spacing>
                default = 0.005r
                </spacing>
                radius    = 0.9r
                thickness = 40p
                fill      = yes
                show_label = yes
                label_font = default
                label_radius = 1.08r
                label_size = 40
                label_parallel = yes
                show_bands = yes
                fill_bands = yes
                </ideogram>""")
        self.ideogram = filename
        return filename

    def write_ticks_config(self, filename='txticks.conf'):
        with open(filename,'w') as f:
            f.write("""show_ticks = yes
                show_tick_labels = yes
                <ticks>
                radius = 1.01r
                color = black
                thickness = 2p
                multiplier = 1e-6
                format = %d
                <tick>
                spacing = 1u
                size = 5p
                </tick>
                <tick>
                spacing = 5u
                size = 10p
                show_label = yes
                label_size = 20p
                label_offset = 10p
                format = %d
                </tick>
                </ticks>""")
        self.ticks = filename
        return filename

    def generate_config(self, ticks = 'txticks.conf', ideogram = 'txideogram.conf', links_and_rules = 'linksAndrules.conf', config='circos.conf'):
        colors = pd.read_table(self.synteny.s_genome.karyotype,header=None,usecols=[2,6],sep=' ').as_matrix()
        self.links_and_rules = links_and_rules
        self.config = config
        self.ideogram,self.ticks = ideogram, ticks
        if hasattr(self, 'ticks'):
            self.write_ticks_config(self.ticks)
        if hasattr(self, 'ideogram'):
            self.write_ideogram_config(self.ideogram)
        with open(self.config,'w') as f:
            f.write("""# circos.conf
                karyotype = %s, %s
                chromosomes_units = 1000000
                chromosomes_display_default = yes
                <<include %s>>
                <<include %s>>
                <<include %s>>
                <image>
                <<include etc/image.conf>>
                </image>
                <<include etc/colors_fonts_patterns.conf>>
                <<include etc/housekeeping.conf>>
                """%(self.synteny.q_genome.karyotype,self.synteny.s_genome.karyotype,self.ideogram,self.ticks,self.links_and_rules))
        with open(self.links_and_rules,'w') as f:
            f.write("""
                <links>
                <link>
                file = %s
                radius = 0.99r
                bezier_radius = 0r
                ribbon = yes
                color = black_a4
                <rules>
                <rule>
                condition = var(intrachr)
                show = no
                </rule>\n"""%(self.synteny.link) + '\n'.join(['<rule>\ncondition = to(%s)\ncolor = %s\n</rule>'%(chrom,color) for chrom,color in map(tuple,colors)]) + '\n</rules>\n</link>\n</links>')

    def run_circos(self, output_dir='./', pdf=False):
        subprocess.call('circos -conf %s -outputfile %s-%s -outputdir %s'%(self.config,self.synteny.q_genome.protID,self.synteny.s_genome.protID,output_dir),shell=True)
        if pdf:
            subprocess.call('convert %s/%s-%s.png %s/%s-%s.pdf'%(os.path.abspath(output_dir),self.synteny.q_genome.protID,self.synteny.s_genome.protID,os.path.abspath(output_dir),self.synteny.q_genome.protID,self.synteny.s_genome.protID),shell=True)

#######################
#### RUN CLI GROUP ####

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'], max_content_width=90)

@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version='0.2')
def circosdrop():
    pass

####################
#### RUN CIRCOS ####

@circosdrop.command()
@click.option('-fi', '--fasta_path', default = './fasta_path/', show_default=True, help='Fasta path containing all of the input genomes. Genome naming must conform to xxx_[protID]_xxx.[fa/fasta].', type=click.Path(exists=False))
@click.option('-gff', '--gff_path', default = './gff_path/', show_default=True, help='Gff path containing all of the gff/gff3 files. Gff naming must conform to: xxx.[protID].[gff/gff3].', type=click.Path(exists=False))
@click.option('-s', '--synteny_path', default = './synteny_path/', show_default=True, help='Path containing synteny files, .unout or .anchors files. *.unout must conform to following pattern: [PAC4GC/PAC2_0].[q_protID]-[PAC4GC/PAC2_0].[s_protID]_5.unout; *.anchors must conform to: [q_protID].[s_protID].[*].anchors. Not neccessary to add files to this path, synteny will be generated if no specification.', type=click.Path(exists=False))
@click.option('-bed', '--bed_path', default = './bed_path/', show_default=True, help='Bed path containing all of the bed files.', type=click.Path(exists=False))
@click.option('-ci', '--circos_inputs', default = './circos_inputs/', show_default=True, help='Path containing all of circos inputs and configuration files.', type=click.Path(exists=False))
@click.option('-co', '--circos_outputs', default = './circos_outputs/', show_default=True, help='Path containing all of circos output images.', type=click.Path(exists=False))
@click.option('-l', '--loci_threshold', default= 4, show_default=True, help='Minimum number of genes in a syntenic block in order to include the block.')
@click.option('-info', '--gene_info', default = 'Name', show_default=True, help='Naming convention for gff file\'s gene name field.', type=click.Choice(['Name', 'gene_name']))
@click.option('-n', '--n_cpus', default = 16, show_default=True, help='Number of cpus used to convert hal 2 maf files.')
def circos_dropper(fasta_path, gff_path, synteny_path, bed_path, circos_inputs, circos_outputs, loci_threshold,gene_info, n_cpus):
    """Visualize many pairwise synteny results. If synteny files are not supplied, conduct pairwise synteny between all included strains."""
    fasta_files = {fasta.split('_')[-2] : fasta for fasta in glob.glob(fasta_path+'/*.fa')+glob.glob(fasta_path+'/*.fasta')}
    gff_files = {gff.split('.')[-2] : gff for gff in glob.glob(gff_path+'/*.gff')+glob.glob(gff_path+'/*.gff3') }
    intersect_keys = set(fasta_files.keys()) & set(gff_files.keys())
    fasta_files = {protID:fasta for protID,fasta in fasta_files.items() if protID in intersect_keys}
    gff_files = {protID:gff for protID,gff in gff_files.items() if protID in intersect_keys}
    genomes = {}
    pairwise_syntenies = []
    print gff_files, fasta_files
    for protID in intersect_keys:
        genomes[protID] = Genome(fasta_files[protID],bed_path+'/'+protID+'.bed',protID,gff_files[protID],gene_info)
        genomes[protID].export_karyotype(circos_inputs+'/'+protID+'.karyotype.txt')
    print genomes
    synteny_files = glob.glob(synteny_path+'/*.unout')+glob.glob(synteny_path+'/*.lifted.anchors')
    synteny_protIDs = []
    if synteny_files:
        for synteny_file in synteny_files:
            if synteny_file.endswith('.unout'):
                coords = reduce(lambda x,y: x+y, sorted([[m.start(0),m.end(0)] for m in re.finditer('PAC2_0|PAC4GC',synteny_file)]))[1::2]
                q_protID, s_prot_ID = map(lambda x: synteny_file[x+1:x+4],coords)
            else:
                q_protID, s_prot_ID = tuple(synteny_file[synteny_file.rfind('/')+1:].split('.')[:2])
            synteny_protIDs.extend([(q_protID, s_prot_ID),(s_prot_ID, q_protID)])
            pairwise_synteny = PairwiseSynteny(genomes[q_protID],genomes[s_prot_ID],synteny_file,loci_threshold=loci_threshold)
            pairwise_synteny.generate_synteny_structure(synteny_path)
            pairwise_syntenies.append(pairwise_synteny)
        remaining_synteny = set(list(combinations(intersect_keys,r=2))) - set(synteny_protIDs)
    else:
        remaining_synteny = list(combinations(intersect_keys,r=2))
    if remaining_synteny:

        def generate_CDS(protID):
            print(protID)
            genomes[protID].extract_CDS()
            return protID
        p = mp.ProcessingPool(n_cpus)
        r = p.amap(generate_CDS,list(set(reduce(lambda x,y: list(x)+list(y),remaining_synteny))))
        r.wait()
        protIDs = r.get()


        def p_synteny(protIDs):
            print(protIDs)
            q_protID, s_prot_ID = protIDs
            pairwise_synteny = PairwiseSynteny(genomes[q_protID],genomes[s_prot_ID],loci_threshold=loci_threshold)
            pairwise_synteny.generate_synteny_structure(synteny_path)
            return pairwise_synteny

        """
        for q_protID, s_prot_ID in combinations(genomes.keys(),r=2):
            pairwise_synteny = PairwiseSynteny(genomes[q_protID],genomes[s_prot_ID],loci_threshold=loci_threshold)
            pairwise_synteny.generate_synteny_structure(synteny_path)
            pairwise_syntenies.append(pairwise_synteny)"""
        r = p.amap(p_synteny,remaining_synteny)#,callback=mycallback) # _async
        r.wait()
        pairwise_syntenies.extend(r.get())
        p.close()


    for pairwise_synteny in pairwise_syntenies:
        pairwise_synteny.synteny_structure_2_link(circos_inputs+'/%s.%s.link.txt'%(pairwise_synteny.q_genome.protID,pairwise_synteny.s_genome.protID))
        circos_obj = Circos(pairwise_synteny)
        circos_obj.generate_config(ticks = circos_inputs+'/txticks.conf', ideogram = circos_inputs+'/txideogram.conf', links_and_rules = circos_inputs+'/linksAndrules.conf', config=circos_inputs+'/circos.conf')
        circos_obj.run_circos(circos_outputs+'/')

@circosdrop.command()
@click.option('-fi', '--fasta_path', default = './fasta_path/', show_default=True, help='Fasta path containing all of the input genomes. Genome naming must conform to xxx_[protID]_xxx.[fa/fasta].', type=click.Path(exists=False))
@click.option('-gff', '--gff_path', default = './gff_path/', show_default=True, help='Gff path containing all of the gff/gff3 files. Gff naming must conform to: xxx.[protID].[gff/gff3].', type=click.Path(exists=False))
@click.option('-bed', '--bed_path', default = './bed_path/', show_default=True, help='Bed path containing all of the bed files.', type=click.Path(exists=False))
@click.option('-w', '--work_dir', default = './', show_default=True, help='Work directory.', type=click.Path(exists=False))
@click.option('-info', '--gene_info', default = 'Name', show_default=True, help='Naming convention for gff file\'s gene name field.', type=click.Choice(['Name', 'gene_name']))
def generate_genomes(fasta_path,gff_path, bed_path, work_dir, gene_info):
    fasta_files = {fasta.split('_')[-2] : os.path.abspath(fasta) for fasta in glob.glob(fasta_path+'/*.fa')+glob.glob(fasta_path+'/*.fasta')}
    gff_files = {gff.split('.')[-2] : os.path.abspath(gff) for gff in glob.glob(gff_path+'/*.gff')+glob.glob(gff_path+'/*.gff3') }
    intersect_keys = set(fasta_files.keys()) & set(gff_files.keys())
    fasta_files = {protID:fasta for protID,fasta in fasta_files.items() if protID in intersect_keys}
    gff_files = {protID:gff for protID,gff in gff_files.items() if protID in intersect_keys}
    genomes = {}
    for protID in intersect_keys:
        genomes[protID] = Genome(fasta_files[protID],bed_path+'/'+protID+'.bed',protID,gff_files[protID],gene_info)
    for genome in genomes:
        pickle.dump(genomes[genome],open('%s%s.p'%(work_dir,protID),'wb'))

@circosdrop.command()
@click.option('-i', '--genome_pickle', help='Input genome pickle.', type=click.Path(exists=False))
def generate_cds(genome_pickle):
    pickle.load(open(genome_pickle,'rb')).extract_CDS()

@circosdrop.command()
@click.option('-i', '--genome_pickle', help='Input genome pickle.', type=click.Path(exists=False))
@click.option('-ci', '--circos_inputs', default = './circos_inputs/', show_default=True, help='Path containing all of circos inputs and configuration files.', type=click.Path(exists=False))
def generate_karyotype(genome_pickle, circos_inputs):
    genome = pickle.load(open(genome_pickle,'rb'))
    genome.export_karyotype(circos_inputs+'/'+genome.protID+'.karyotype.txt')

@circosdrop.command()
@click.option('-w', '--work_dir', default = './', show_default=True, help='Work directory.', type=click.Path(exists=False))
@click.option('-g','--genomes',help='Genome pickles, space delimited.',multiple=True)
def pair_genomes(work_dir,genomes):
    genomes = {genome.protID:genome for genome in pickle.load(open(genomes,'rb'))}
    pairs = set(list(combinations(genomes.keys(), r=2)))
    for pair in pairs:
        pickle.dump((genomes[pair[0]],genomes[pair[1]]),open('%s%s.%s.p' %(work_dir,pair[0],pair[1]),'wb'))

@circosdrop.command()
@click.option('-g', '--genome_pair', help='Pair of input genome pickles.', type=click.Path(exists=False))
@click.option('-l', '--loci_threshold', default= 4, show_default=True, help='Minimum number of genes in a syntenic block in order to include the block.')
@click.option('-s', '--synteny_path', default = './synteny_path/', show_default=True, help='Path containing synteny files, .unout or .anchors files. *.unout must conform to following pattern: [PAC4GC/PAC2_0].[q_protID]-[PAC4GC/PAC2_0].[s_protID]_5.unout; *.anchors must conform to: [q_protID].[s_protID].[*].anchors. Not neccessary to add files to this path, synteny will be generated if no specification.', type=click.Path(exists=False))
@click.option('-w', '--work_dir', default = './', show_default=True, help='Work directory.', type=click.Path(exists=False))
def generate_synteny(genome_pair, loci_threshold, synteny_path, work_dir):
    genome1, genome2 = pickle.load(open(genome_pair,'rb'))
    pairwise_synteny = PairwiseSynteny(genome1, genome2, loci_threshold=loci_threshold)
    pairwise_synteny.generate_synteny_structure(synteny_path)
    pickle.dump(pairwise_synteny, open('%s%s.%s.p' % (work_dir, pairwise_synteny.q_genome.protID, pairwise_synteny.s_genome.protID), 'wb'))

@circosdrop.command()
@click.option('-p', '--pairwise_pickle', help='Input genome pickle 1.', type=click.Path(exists=False))
@click.option('-ci', '--circos_inputs', default = './circos_inputs/', show_default=True, help='Path containing all of circos inputs and configuration files.', type=click.Path(exists=False))
@click.option('-co', '--circos_outputs', default = './circos_outputs/', show_default=True, help='Path containing all of circos output images.', type=click.Path(exists=False))
def display_circos(pairwise_pickle,circos_inputs, circos_outputs):
    pairwise_synteny = pickle.load(open(pairwise_pickle,'rb'))
    pairwise_synteny.synteny_structure_2_link(circos_inputs + '/%s.%s.link.txt' % (pairwise_synteny.q_genome.protID, pairwise_synteny.s_genome.protID))
    circos_obj = Circos(pairwise_synteny)
    circos_obj.generate_config(ticks=circos_inputs + '/txticks.conf', ideogram=circos_inputs + '/txideogram.conf',
                               links_and_rules=circos_inputs + '/linksAndrules.conf',
                               config=circos_inputs + '/circos.conf')
    circos_obj.run_circos(circos_outputs + '/')



#### RUN CLI ####

if __name__ == '__main__':
    circosdrop()
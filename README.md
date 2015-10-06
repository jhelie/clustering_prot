# clustering_prot
Python script to analyse clustering of proteins.
To print the help below: ```python clustering_prot.py --help```

```
[ DESCRITPION ]

This script calculates the evolution of proteins clustering status.
	
It produces 3 outputs ((a) and (b) require --groups to be specified, see note 5):
 (a) 2D plots: time evolution of the cluster size each protein is involved in
 (b) 1D plots: time evolution of the % of protein represented by each group
 (c) stability statistics: max nb of consecutive frames each group existed for [TO DO]

The metrics can be calculated for a single frame or for an entire trajectory - and in case
a trajectory is supplied the data for individual frame snapshots can also be produced at a
frequency specified by the option -w (if you don't specify an argument to -w snapshots
will only be written for the first and last frame). Note that writing files to disk will
considerably slow down the execution of the script.

Detection of transmembrane protein clusters
-------------------------------------------
Two clustering algorithms can be used to identify protein clusters.
->Connectivity based (relies on networkX module):
  A protein is considered in a cluster if it is within a distance less than --nx_cutoff
  from another protein. This means that a single protein can act as a connector between
  two otherwise disconnected protein clusters.
  This algorithm can be ran using either the minimum distante between proteins (default, 
  --algorithm 'min') or the distance between their center of geometry (--algorithm 'cog').
  The 'min' option scales as the square of the number of proteins and can thus be very
  slow for large systems.

->Density based (relies on the sklearn module and its implementation of DBSCAN):
  A protein is considered in a cluster if is surrounded by at least --db_neighbours other
  proteins within a radius of --db_radius.
  This density based approach is usually less suited to the detection of protein
  clusters but as a general rule the more compact the clusters, the smaller --db_radius
  the higher --db_neighbours can be - for details on this algorithm see its online
  documentation.
  This algorithm is selected by setting the --algorithm option to 'density'.

The identified protein clusters are considered to be transmembrane only if the closest
lipid headgroup neighbours to the cluster particles are all within the same leaflet.
In addition to the sizes identified, size groups can be defined - see note 7.

VMD visualisation
-----------------
The perturbation calculated can be visualised with VMD either for a single frame or an
entire trajectory. Note that in the case of a trajectory only the processed frame will be
annotated (every 10 by defaults) - you might want to pre-process your trajectory to remove
frames and then set the -t option to 1 when running the script.
 ->frame:
   The thickness or order parameter info is stored in the beta factor column. Just open
   the PDB file with VMD and choose Draw Style > Coloring Method > Beta.

 ->trajectory:
   The thickness or order parameter info is stored in a .txt file in folders '/3_VMD/'.
   To load it into your trajectory in VMD use the appropriate routine of the script
   script 'vmd_parser.tcl' (see https://github.com/jhelie/vmd_parser).


[ REQUIREMENTS ]

The following python modules are needed :
 - MDAnalysis
 - matplotlib
 - networkX (if --algorithm set to 'min' or 'cog')
 - sklearn (if --algorithm set to 'density')


[ NOTES ]

1. It's a good idea to trjconv the xtc first and only outputs the proteins,  as the 
   script will run MUCH faster. Also, use the -pbc mol option.	

2. The script tracks the accummulation of proteins on each leaflet, but detailed clustering
   information is only calculated for transmembrane proteins.
   
   Identification of the bilayer leaflets can be further controlled via 3 ways:
   (a) beads
    By default, the particles taken into account to define leaflet are:e
    -> name PO4 or name PO3 or name B1A
   
    Note that only lipids which contain one of the beads mentioned in the selection string
    will be taken into account. If you wish to specify your own selection string (e.g. to
    choose different beads or add a bead not in the default list in order to take into
    account a particular lipid specie) you can do so by supplying a file via the --beads
    option. This file should contain a single line that can be passed as the argument
    to MDAnalysis selectAtoms() routine and should not contain any quotation marks, e.g.:
     -> name PO4 or name PO3 or name B1A or name AM1
        
   (b) leaflet finding method
    By default leaflets are identified using the MDAnalysis LeafletFinder routine and the
    the optimum cutoff to identify 2 lipids groups is determined using the optimize_cutoff
    routine.
    This optimisation process can take time in large systems and you can specify your own
    cutoff value to skip this step. For instance to use the default 15 Angstrom cutoff
    directly (without optimising):
     -> '--leaflets 15'
   
    In very large systems (more then ~50,000 phospholipids) LeafletFinder (or rather the
    networkX module that it relies on) can fail. To  avoid this you can choose not to use
    this routine by specifying:
     -> '--leaflets large'
    In this case lipids whose headgroups z value is above the average lipids z value will
    be considered to make up the upper leaflet and those whose headgroups z value is below
    the average will be considered to be in the lower leaflet.
    This means that the bilayer should be as flat as possible in the 1st frame of the xtc
    file supplied in order to get a meaningful outcome. 

	NOTE: By default the gro file is only used as a topology file and the 1st frame of the
	xtc is used to identify leaflets. If you wish to use the gro file instead, for instance
	in the case that the 1st frame of the xtc is not flat, you need to specify the --use_gro
	flag: be warned that this might take a few minutes longer on large systems.

   (c) flipflopping lipids
    In case lipids flipflop during the trajectory, a file listing them can be supplied
    with the --flipflops option. Each line of this file should follow the format:
     -> 'resname,resid,starting_leaflet,z_bead'
    where starting_leaflet is either 'upper' or 'lower' - e.g. 'POPC,145,lower,PO4'. The
    z_bead is used to track the position of the lipid.
    If flipflopping lipids are not specified they may add significant noise to results.   

3. Proteins are detected automatically but you can specify an input file to define your
   own selection with the -p option.
   In this case the supplied file should contain on each line a protein selection string
   that can be passed as the argument of the MDAnalysis selectAtoms() routine - for 
   instance 'bynum 1:344'.

4. Clusters statistics can be binned into size groups. The size groups are defined by
   supplying a file with the --groups option, whose lines all follow the format:
    -> 'lower_size,upper_size, colour'

   Size groups definition should follow the following rules:
    -to specify an open ended group use 'max', e.g. '3,max,colour'
    -groups should be ordered by increasing size and their boundaries should not overlap
    -boundaries are inclusive so you can specify one size groups with 'size,size,colour'
    -colours must be specified for each group (see (c) in note 5)
    -any cluster size not falling within the specified size groups will be labeled as
     'other' and coloured in grey (#E6E6E6)
	
5. The colours associated to each TM cluster size identified are based on the matplotlib
   'jet' colour map by default. You can specify your own colours as follows:
   (a) Colour of individual cluster sizes
    Colours of individual cluster sizes use the matplotlib 'jet' colour scheme and cannot 
    be modified. 

   (b) Colour of size groups
	See note 4 above: exact colour control can be achieved for each size by using size
	groups of unitary size, so if you want to control individual cluster sizes colours just
	specify the relevant group file.

   (c) Colour definition
    Colours can be specified using single letter code (rgbcmykw), hex code  or the name of
    a colour map (see the matplotlib website for a list of the available colour maps).
    In case a colour map is used, its name must be specified as the colour for each lipid
    specie or size group.

6. The script automatically detects proteins and, by default, groups them into species labelled
   'A', 'B', etc... based on sequences. Oligomers, ie proteins made up of the repeat of a
   same sequence, are not detected by default and just considered as a protein with the overall
   sequence which can effectively decrease your sampling of contacts on heatmaps representing
   residues interactions.
   The --species flag allows to specify flag to name proteins and address the oligomer issue.
   Each line of this file should follow the format (without the quotes):
    ->'name,multiplicity,sequence'
    
   with the multiplicity referring to the number of times the sequence appears. So, in the
   case of, say, a trimer 'multiplicity' can be set to 3 with 'sequence' corresponding to
   the sequence of a monomer. For proteins which do not involve the repetition of a base
   unit 'multiplicity' should be 1 and 'sequence' the complete protein sequence.
   1 letter code should be used for sequence.

7  The --nx_cutoff option is mandatory and allows to control the distance for which proteins are
   considered to be in contact.
 
   In case the 'min' option for the --algorithm flag is specified, the --nx_cutoff flag should be
   set to a float value (typically 6 or 8 Angstrom in CG).
   
   In case the 'cog' option for --algorithm is used (which should always be the case for rigid TM
   proteins as it is MUCH faster) the --nx_cutoff can either be set to a float value (Angstrom) or
   used to specify a text file. This latter option is particularly useful if the system contains
   several protein species of very different sizes. In this case a cutoff specific to each type of
   interactions must be supplied in the txt file, with each line following the format:
    -> 'specie1_index,specie2_index,cutoff_distance_between_cog1_and_cog2'
    
    where 'specie1_index' and 'specie2_index' correspond to 0-based indices referring to a specie of
    proteins present in the system (numbered according to their appearance in the pdb or gro file).
    For instance, for a system containing two types of proteins, the following would be acceptable:
     '0,0,60'
     '1,1,100'
     '0,1,80'
   
	NB: this means that the systems MUST be built so that the different protein species are listed
	sequentially in the gro/pdb files (they cannot alternate) ant that you know the order in which
	the species are listed. Doing otherwise would generally be poor practice and lead to complex
	topology files anyway.
 

[ USAGE ]

Option	      Default  	Description                    
-----------------------------------------------------
-f			: structure file [.gro] (required)
-x			: trajectory file [.xtc]
-o			: name of output folder
-b			: beginning time (ns)
-e			: ending time (ns)	
-t 		10	: process every t-frames
-w			: write snapshots every [w] processed frames (see 'DESCRIPTION')

Lipids identification (see note 2)
-----------------------------------------------------
--beads			: leaflet identification technique, see note 2(a)
--flipflops		: input file with flipflopping lipids, see note 2(c)
--leaflets	optimise: leaflet identification 'optimise', 'large', 'no' or float, see note 2(b)
--use_gro		: use gro file instead of xtc, see note 2(b)

Proteins properties
-----------------------------------------------------
--species		: file defining name,multiplicity and sequence of proteins, see note 6
--groups		: cluster groups definition file, see note 4 [BETA]
--res_contact	8	: cutoff to consider contacts between residues c.o.g (Angstrom)
--res_show	1	: show interactions for residues for at least that many % of contacts between proteins
--algorithm	cog	: 'cog','min' or 'density', see 'DESCRIPTION'
--nx_cutoff 		: networkX cutoff distance for protein-protein contact (Angstrom or file), see note 7
--db_radius 	20	: DBSCAN search radius (Angstrom)
--db_neighbours	3	: DBSCAN minimum number of neighbours within a circle of radius --db_radius	
 
Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
```

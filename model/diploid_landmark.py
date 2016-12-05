'''Nuclear Landmark model for budding yeast'''
import os
import sys
import numpy as np
import datetime
import string
import math
import IMP
import IMP.core
import IMP.atom
import IMP.display
import IMP.algebra
import unittest
from StringIO import StringIO
import time
import argparse
import csv

#parse arguments
parser = argparse.ArgumentParser(description='Run volume exclusion model')
parser.add_argument('seed', type=int, help='seed')
parser.add_argument('name', type=str, help='model name')
parser.add_argument('chrlen', type=str, help='chromosome lengths')
parser.add_argument('chrcen', type=str, help='centromeres')
parser.add_argument('size', type=float, help='nuclear radius')
args = parser.parse_args()

# Chromosome lengths and list of chromosomes
chrlens = list(csv.reader(open(args.chrlen,'rb'),delimiter='\t'))
chr_seq = {}
chain_list = []
for i in range(0,len(chrlens)):
    chr_seq[chrlens[i][0]] = int(chrlens[i][1])
    chain_list.append(chrlens[i][0])

# Chromosome centromeres
chrcens = list(csv.reader(open(args.chrcen,'rb'),delimiter='\t'))
chr_cen = {}
for i in range(0,len(chrcens)):
    chr_cen[chrcens[i][0]] = (int(chrcens[i][1]) + int(chrcens[i][2]))/2

# Chromosome names for PDB
chr_pdb = {}
for i in range(0,len(chrlens)):
    if i < 16:
        chr_pdb[chrlens[i][0]] = "c" + str(i+1).zfill(2) + " " + string.ascii_uppercase[i]
    else:
        chr_pdb[chrlens[i][0]] = "c" + str(i+1).zfill(2) + " " + string.ascii_lowercase[i-16]

# parameters
sep = 3200       # 3200 bp separation
nuclear_rad = 1000.0 * args.size
nucleolus_pos = -1200.0 * args.size
nucleolus_in_rad = 750.0 * args.size
nucleolus_ex_rad = 1600.0 * args.size
centro_pos = -750.0 * args.size
centro_rad = 300.0 * args.size
envelope_thick = 50 * args.size

# Set the seed and the random state
random_state = np.random.RandomState(seed=args.seed)

basedir = os.path.dirname(os.path.realpath(__file__))
outname = os.path.join(basedir, args.name, str(args.seed))
print 'nucleus R={}, cen {}, R={}, nucleolus {}, Rin={}, Rex={}, tel {} from NE, seed={}'.format(nuclear_rad,centro_pos,centro_rad,nucleolus_pos,nucleolus_in_rad,nucleolus_ex_rad,envelope_thick,args.seed)
print 'output pdb:', outname
if os.path.exists(outname):
    sys.exit(0)

#------------------------------------------------------------

t1=time.time()
chr_bead = {}    # number of beads for each chromosome
nbead = 0
bead_start = {}  # bead label starts of a chr
for i in chr_seq.keys():
    n = chr_seq[i]/sep + 1
    chr_bead[i] = n
    nbead = nbead + n
    bead_start[i] = nbead - n

#print bead_start
rdnaStart = bead_start[chain_list[11]] + 140 #begin rDNA
rdnaEnd = bead_start[chain_list[11]] + 147 #not included as rDNA
rdna1 = rdnaStart + 2 #last bead 1st chain
rdna2 = rdnaStart + 3 #begin 2nd chain
rdnaStart2 = bead_start[chain_list[27]] + 140 #begin rDNA
rdnaEnd2 = bead_start[chain_list[27]] + 147 #not included as rDNA
rdna3 = rdnaStart2 + 2 #last bead 1st chain
rdna4 = rdnaStart2 + 3 #begin 2nd chain

startChr_pdb = [] #indexing starts from 0
n = 0
startChr_pdb.append(n)
for i in chain_list:
    n += chr_bead[i]
    startChr_pdb.append(n)

#---------------------------------------------------------

def bead_id(chr,gpos):
    '''Given chromosome id and genome position, returns bead id'''
    for i in range(chr_bead[chr]):
        if gpos >= i*sep and gpos < (i+1)*sep:
            beadnum = i + bead_start[chr]
            break
    return beadnum

def find_chromosome(bid):
    """ Returns a chromosome id given a bead number"""
    for i in chr_seq.keys():
        if bid < bead_start[i] + chr_bead[i]\
        and bid >= bead_start[i]:
            chrid=i
            break
    return chrid

def find_bead_in_chr(bid):
    """ Returns a chromosome and bead_order and mid genome position given a beadnum"""
    for i in chr_seq.keys():
        if bid < bead_start[i] + chr_bead[i]\
        and bid >= bead_start[i]:
            order = bid - bead_start[i] + 1 #order starts from 1
            genpos = order*sep - sep/2
            break
    return i,order,genpos

def pdboutput(name):
    pdb=[[] for i in range(nbead)]
    #-----------------------------
    cen_pos=[]
    for k in chr_seq.keys():
        j= bead_id(k,chr_cen[k])
        cen_pos.append(j)
    #------------------------------------
    for i in range(nbead):
        p0=IMP.core.XYZR(chain.get_particle(i))
        chr=find_chromosome(i)
        #pdb[i].append('ATOM') 
        if i in cen_pos:       
            pdb[i].append(' CEN')  #1
        elif rdnaStart<=i<=rdnaEnd or rdnaStart2<=i<=rdnaEnd2:
            pdb[i].append('rDNA') 
        elif i == bead_start[chr]:
            pdb[i].append(' L  ') 
        elif i == (bead_start[chr]+chr_bead[chr]-1):
            pdb[i].append(' R  ') 
        else:
            pdb[i].append(' O  ')
        #pdb[i].append('L')
        chr_num=filter(lambda k: bead_start[k]<=i, chr_seq.keys())[-1]
        pdb[i].append(chr_pdb[chr_num])  #2
        pdb[i].append(i-bead_start[chr_num]+1) #3
        pdb[i].append(p0.get_x()) #4
        pdb[i].append(p0.get_y()) #5
        pdb[i].append(p0.get_z()) #6
        #pdb[i].append(15)
        pdb[i].append(chr_num)  #7
    ###sort the file by chromosome order
    sorted_pdb=[]
    for i in chain_list:
        for j in range(len(pdb)):
            if pdb[j][6]==i:
                sorted_pdb.append(pdb[j])
    #insert sorted number
    for i in range(len(sorted_pdb)):
        sorted_pdb[i].insert(0,i+1)

    name=str(name)+'.pdb'
    #------------------------------------------------
    out=open(name,'w')
    for l in sorted_pdb:
        out.write("ATOM %6i %4s %5s %3i     %7.1f %7.1f %7.1f %s\n"\
        %(l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7]))
    out.close()

def mdstep(t,step):
    o = IMP.atom.MolecularDynamics()
    o.set_model(m)
    md = IMP.atom.VelocityScalingOptimizerState(xyzr,t,10)  # replace 300 K with 500 K
    o.add_optimizer_state(md)
    #print 'optimizing with temperature',t,'and',step,'steps'
    s=o.optimize(step)
    o.remove_optimizer_state(md)
    #print 'MD',step,'steps done @',datetime.datetime.now()
    return s

def cgstep(step):
    o = IMP.core.ConjugateGradients()
    o.set_model(m)
    f=o.optimize(step)
    #print 'CG',step,'steps done @',datetime.datetime.now()
    return f

#___________________________ IMP starts _____________________________________
#IMP.set_check_level(IMP.NONE)
IMP.set_log_level(IMP.SILENT)
m = IMP.Model()
r = 15.0
lb = 30.0 # length of bond
kbend=0.2 
contact_dict = {}
xyzr = IMP.core.create_xyzr_particles(m,nbead,r)
chain = IMP.container.ListSingletonContainer(xyzr)
# First beads
corner1=IMP.algebra.Vector3D(-nuclear_rad,-nuclear_rad,-nuclear_rad)
corner2=IMP.algebra.Vector3D(nuclear_rad,nuclear_rad,nuclear_rad)
box=IMP.algebra.BoundingBox3D(corner1,corner2)
rdummy=int(random_state.rand() * 10000)
for i in range(rdummy):
    ranvec = IMP.algebra.get_random_vector_in(box)
#----------------------------------------------------------  
#print nbead
for i in range(nbead):
    p0 = chain.get_particle(i)
    IMP.atom.Mass.setup_particle(p0,1)
    p = IMP.core.XYZR(p0)
    coor = IMP.algebra.get_random_vector_in(box)
    p.set_coordinates(coor)
    #ch,b,gpos = find_bead_in_chr(i)
    #print ch,b,pdbOrder(ch,gpos)+1,i

#sys.exit()
#---------------------------------------------------------------------------------
#print 'Setting up restraints'
# Create bonds for consecutive beads in a string
bonds = IMP.container.ListSingletonContainer(m)
for id in chr_seq.keys():
    istart = bead_start[id]
    iend = istart + chr_bead[id]
    IMP.atom.Bonded.setup_particle(chain.get_particle(istart))
    for i in range(istart + 1,iend):
        if i != rdna2 and i != rdna4:
            bp = IMP.atom.Bonded.decorate_particle(chain.get_particle(i-1))
            bpr = IMP.atom.Bonded.setup_particle(chain.get_particle(i))
            b = IMP.atom.create_custom_bond(bp, bpr, lb, 2)
            bonds.add_particle(b.get_particle())
        elif i == rdna4:
            IMP.atom.Bonded.setup_particle(chain.get_particle(rdna4))
        else:
            IMP.atom.Bonded.setup_particle(chain.get_particle(rdna2))

# Restraint for bonds
bss = IMP.atom.BondSingletonScore(IMP.core.Harmonic(0,1))
br = IMP.container.SingletonsRestraint(bss, bonds)
m.add_restraint(br) #0

# Set up excluded volume
evr = IMP.core.ExcludedVolumeRestraint(chain)
m.add_restraint(evr) #1


# Set up cap
center = IMP.algebra.Vector3D(0,0,0)
center_nucleolus = IMP.algebra.Vector3D(nucleolus_pos,0,0)
not_rDNA = IMP.container.ListSingletonContainer(m)
for i in range(nbead):
    if not ((i >= rdnaStart and i < rdnaEnd) or (i >= rdnaStart2 and i < rdnaEnd2)):
        p = chain.get_particle(i)
        not_rDNA.add_particle(p)
ubcell = IMP.core.HarmonicUpperBound(nuclear_rad,1.0)
sscell = IMP.core.DistanceToSingletonScore(ubcell,center)
rcell = IMP.container.SingletonsRestraint(sscell,chain)
m.add_restraint(rcell) #2

# centromeres in radius 300 @-700
centro = IMP.algebra.Vector3D(centro_pos,0,0)
listcentro = IMP.container.ListSingletonContainer(m)
for k in chr_seq.keys():
    j = bead_id(k,chr_cen[k])
    pcen = chain.get_particle(j)
    listcentro.add_particle(pcen)

ubcen = IMP.core.HarmonicUpperBound(centro_rad,1.0)
sscen = IMP.core.DistanceToSingletonScore(ubcen,centro)
rcentro = IMP.container.SingletonsRestraint(sscen,listcentro)
m.add_restraint(rcentro) #4

print 'High temp MD..'
mdstep(1000000,500)
mdstep(500000,500)
mdstep(300000,500)
mdstep(100000,500)
mdstep(5000,500)
score=cgstep(1000)
print 'before telo: ',score

# Telomeres near nuclear envelope thickness 50
telo = IMP.container.ListSingletonContainer(m)
#galBead =  bead_id(galpos[0],galpos[1])
#telo.add_particle(chain.get_particle(galBead))
for k in chr_seq.keys():
    j1 = bead_start[k]
    pt = chain.get_particle(j1)
    telo.add_particle(pt)
    j2 = j1 - 1 + chr_bead[k]
    pt = chain.get_particle(j2)
    telo.add_particle(pt)
tlb = IMP.core.HarmonicLowerBound(nuclear_rad-envelope_thick,1.0)
sst = IMP.core.DistanceToSingletonScore(tlb,center)
rt = IMP.container.SingletonsRestraint(sst,telo)
m.add_restraint(rt) #5

# outside centro sphere
#lbcen = IMP.core.HarmonicLowerBound(centro_rad,1.0)
#sstc = IMP.core.DistanceToSingletonScore(lbcen,centro)
#rtc = IMP.container.SingletonsRestraint(sstc,telo)
#m.add_restraint(rtc) #6

# rDNA chr12 near nucleolus
rDNA = IMP.container.ListSingletonContainer(m)
rDNA.add_particle(chain.get_particle(rdna1))
rDNA.add_particle(chain.get_particle(rdna3))
ub_bn2 = IMP.core.HarmonicLowerBound(nucleolus_in_rad,0.5)  #can also try lowerbound nucleolus_rad-rncutoff
mindts2 = IMP.core.DistanceToSingletonScore(ub_bn2,center_nucleolus)
nucleolir2 = IMP.container.SingletonsRestraint(mindts2,rDNA)
m.add_restraint(nucleolir2) #7

rDNA2 = IMP.container.ListSingletonContainer(m)
rDNA2.add_particle(chain.get_particle(rdna2))
rDNA2.add_particle(chain.get_particle(rdna4))
ub_bn3 = IMP.core.HarmonicLowerBound(nucleolus_ex_rad,0.5)  #can also try lowerbound nucleolus_rad-rncutoff
mindts3 = IMP.core.DistanceToSingletonScore(ub_bn3,center_nucleolus)
nucleolir3 = IMP.container.SingletonsRestraint(mindts3,rDNA2)
m.add_restraint(nucleolir3) #7

# Set up Nucleolis #
lbn0 = IMP.core.HarmonicUpperBound(nucleolus_ex_rad,0.5)
ssn0 = IMP.core.DistanceToSingletonScore(lbn0,center_nucleolus)
rn0 = IMP.container.SingletonsRestraint(ssn0,not_rDNA)
m.add_restraint(rn0) #3
#-------------------

print 'High temp MD in nuc ...'
mdstep(500000,5000)
mdstep(300000,5000)
mdstep(5000,10000)
score=cgstep(500)
print 'before angle',score

# Angle Restraint
angle = math.pi
angle_set=[]
noangle=[i for i in bead_start.values()]  #do not apply angle restraints
noangle.append(rdna1)
noangle.append(rdna2)
noangle.append(rdna3)
noangle.append(rdna4)

for i in range(nbead-1):
    ieval = i+1
    if ieval in noangle:
        continue
    elif i in noangle:
        continue
    else:
        d1 = chain.get_particle(i-1)
        d2 = chain.get_particle(i)
        d3 = chain.get_particle(i+1)
        pot = IMP.core.Harmonic(angle,kbend)
        ar = IMP.core.AngleRestraint(pot,d1,d2,d3)
        m.add_restraint(ar)
        angle_set.append(ar)

mdstep(50000,500)
mdstep(25000,500)
mdstep(20000,1000)
mdstep(10000,1000)
mdstep(5000,3000)
mdstep(2000,5000)
mdstep(1000,7000)
mdstep(500,10000)
score=cgstep(2500)

print 'angle: %.1f '%(score)
#-----------------------
for i in angle_set:
    m.remove_restraint(i)
score=cgstep(1000)
print 'Final score:%.1f' %(score)

pdboutput(outname)

#-------------------------

#mdstep(1000,1000)
#score=cgstep(1000)
#print '\nFinal score with remove angle restraint is: ',score

#name='final_without_angle'
#output(chain,nbead,name)
t2=time.time()

print 'time spend is ', t2-t1, ' s'

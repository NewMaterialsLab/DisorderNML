#!/usr/bin/env python 
# DisorderNML is an application to explore the strcutural, electronic and optical properties 
# of off-stoichiometric materials combined with SOD (https://github.com/gcmt-group/sod) and VASP outputs. 
################################################################################
#      Copyright Jose J. Plata, Antonio M. Marquez and Javier Fdez Sanz        #     
#                               (2020)                                         #
#                                                                              #
# This is free software: you can                                               #
# redistribute it and/or modify it under the terms of the GNU General Public   #
# License as published by the Free Software Foundation, either version 3 of    #
# the License, or (at your option) any later version.                          #
# This program is distributed in the hope that it will be useful, but WITHOUT  #
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        #
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for    #
# more details.                                                                #
# You should have received a copy of the GNU General Public License along with #
# this program.  If not, see <http://www.gnu.org/licenses/>.                   #
#                                                                              #
################################################################################

import os
import sys
import re
import numpy as np 
import xml 
import xml.etree.ElementTree as ET  
import math 
import cmath 
import argparse
from operator import itemgetter
from scipy.integrate import simps
from scipy.integrate import trapz
from scipy.integrate import cumtrapz


KB = 0.000086173303 

def getCalcbins(distlist,prob,binw,rdfm):

	Nbins = rdfm/binw	
	freq = [0.0]*(int(Nbins)+1)
	for dist in distlist:
		if dist < rdfm:
			ind = int((dist // binw) + 1)
			freq[ind] = freq[ind]+(1.0*prob/binw)
	return freq	
		 
def getRdfList(latt,elements,struct,rdfA,rdfB,dim):
	#Calculating indexes in the structure
	#Index A
	count = 0 
	for i in range(rdfA):
		count = count + int(elements[i])
	inA1 = count
	inA2 = count+int(elements[rdfA])
	#Index B
	count = 0 
	for i in range(rdfB):
		count = count + int(elements[i])
	inB1 = count
	inB2 = count+int(elements[rdfB])

	distlist = []
	vec = np.zeros(3)
	vec2 = np.zeros(3)
	for i in range(inA1,inA2):
		#Centering the cell in atom A
		centStruc = centering(struct,i)
		for j in range(inB1,inB2):
			if i < j:
				#Calculate distances in al supercells. Loop over neighbours
				translation = np.zeros(3)
				for o in range(-dim[0],dim[0]+1):
					translation[0] = float(o) 
					for p in range(-dim[0],dim[0]+1):
						translation[1] = float(p) 
						for q in range(-dim[0],dim[0]+1):
							translation[2] = float(q) 
							vec = centStruc[j]+translation
							#Convert Direct in Cartesian
							vec2 = np.matmul(vec,latt)
							distlist.append(np.linalg.norm(vec2))
	distlist.sort()
	return distlist

def centering(struct,ind):
	coord = np.zeros((len(struct),3))
	for i in range(len(struct)):
		coord[i] = struct[i]-struct[ind]
	return coord	

def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		return False

def abcCalc(latt):
	dist = [0.0]*3
	for i in range(3):
		dist[i] = np.linalg.norm(latt[i])
	return dist[0], dist[1], dist[2]

def readPoscar(name):
	os.chdir(name+"/DOS")
	if not os.path.isfile("POSCAR"):
		print("File path {} does not exist. Exiting...".format(name))
		sys.exit()
	f = open("POSCAR", "r")
	lines = f.readlines()
	f.close()
	#Read factor
	chain=lines[1].split()
	factor=float(chain[0])
	#Read Lattice
	latt = np.zeros((3,3))
	for i in range(2,5):
		chain=lines[i].split()
		x = np.array(chain)
		latt[i-2] = x.astype(np.float)
	#Facro*Lattice
	latt = factor * latt
	# Check POSCAR version <4 or >4? If line 6 has floats  vPOSCAR < 4 
	sum_index = 5
	s_index = 7
	chain=lines[5].split()
	v4 = is_number(chain[0])
	if not v4:
		s_index = s_index+1
		sum_index = sum_index+1
    # Selective Dynamics?
	chain=lines[s_index].split()
	NoSelDyn = is_number(chain[0])
	if not NoSelDyn:
		s_index = s_index+1
	#Total atoms summation  
	elements=lines[sum_index].split()
	totalAtom=0
	for at in elements:
		totalAtom = totalAtom + int(at)
	coord = np.zeros((totalAtom,3))
	for i in range(s_index,s_index+totalAtom):
		chain=lines[i].split()
		x = np.array(chain)
		coord[i-s_index] = x[0:3].astype(np.float)
	os.chdir("../..")
	return latt,elements,coord

def gnuHist():
	f= open("hist.plt","w+")
	f.write("""
#!/usr/bin/gnuplot
#
# AUTHOR: Jose J. Plata
# gnuplot 4.6 patchlevel 6

reset

# eps
set terminal postscript eps size 8.6cm,6.437cm  enhanced color \
   font 'Helvetica,18' linewidth 1
#set output 'scell.eps'
set output '| epstopdf --filter --outfile=energy.pdf'


#line styles
set style line 1 lt 1 lc rgb '#0072bd' # light red

# Axes
set style line 101 lc rgb '#000000' lt 1
set tics  in scale 0.75
set mxtics 10
set xtics 1
show mxtics
set mytics 2
show mytics


set ylabel 'Counts'
set xlabel '{{/Symbol D}E (eV)}'

#set xrange [-1:6]
#set yrange [0:350]

unset key
binwidth = 0.1
binstart = 0
load 'hist.fct'

plot 'EnergyHistogram.dat' i 0 @hist ls 1
	""")

	f.close()
	f= open("hist.fct","w+")
	f.write("""
# hist.fct
# gnuplot macro for providing a functionality similar to the hist() function in octave
# Note, that the variables binwidth and binstart has to be set before calling this function
# AUTHOR: Hagen Wierstorf

# set width of single bins in histogram
set boxwidth 0.9*binwidth
# set fill style of bins
set style fill solid 0.5
# define macro for plotting the histogram
hist = 'u (binwidth*(floor(($1-binstart)/binwidth)+0.5)+binstart):(1.0) smooth freq w boxes'
	""")
	f.close()

def gnuSpctra():

	f= open("spectra.plt","w+")
	f.write("""
#!/usr/bin/gnuplot
#
# AUTHOR: Jose J. Plata
# gnuplot 4.6 patchlevel 6

reset

# Preparing terminal, fonts, plot size and output format
set terminal postscript eps size 8.6cm,6.437cm  enhanced color \
   font 'Helvetica,18' linewidth 1
set output '| epstopdf --filter --outfile=2layer_dos.pdf'


#Line styles
set style line 1 lc rgb '#000000' lt 2  # --- black
set style line 2 lc rgb '#dd181f' lt 1 lw 2 pt 7  # --- red
set style line 11 lt 1 lc rgb '#0072bd' lw 2 pt  7  ps 1.5 # blue
set style line 12 lt 1 lc rgb '#77ac30' lw 2 pt  7  ps 1.5# orange

# Axes and tics
set style line 101 lc rgb '#000000' lt 1
set tics  in scale 0.75
set mxtics 2
set xtics 1
show mxtics
set mytics 2
show mytics


set ylabel '{{/Symbol e}_{2}({/Symbol w})}' 
set xlabel '{E (eV)}'

set xrange [0:6]

plot "gsSpectra.dat"          u 1:2       w l ls 11 ti gs,\
	 "avSpectra.dat"          u 1:2       w l ls 12 ti av
	""")
	f.close()

def gnuDOS():

	f= open("dos.plt","w+")
	f.write("""
#!/usr/bin/gnuplot
#
# AUTHOR: Jose J. Plata
# gnuplot 4.6 patchlevel 6

reset

# Preparing terminal, fonts, plot size and output format
set terminal postscript eps size 8.6cm,6.437cm  enhanced color \
   font 'Helvetica,18' linewidth 1
set output '| epstopdf --filter --outfile=2layer_dos.pdf'


#Line styles
set style line 1 lc rgb '#000000' lt 2  # --- black
set style line 2 lc rgb '#dd181f' lt 1 lw 2 pt 7  # --- red
set style line 11 lt 1 lc rgb '#0072bd' lw 2 pt  7  ps 1.5 # blue
set style line 12 lt 1 lc rgb '#77ac30' lw 2 pt  7  ps 1.5# orange

# Axes and tics
set style line 101 lc rgb '#000000' lt 1
set tics  in scale 0.75
set mxtics 2
set xtics 1
show mxtics
set mytics 2
show mytics


set ylabel '{{/Symbol e}_{2}({/Symbol w})}' 
set xlabel '{E (eV)}'

set xrange [0:6]

plot "gsDOS1spin.dat"          u 1:2       w l ls 11 ti gs,\
	 "avDOS1spinXXXK.dat"          u 1:2       w l ls 12 ti av
	""")
	f.close()

def getList(filepath):
	f = open(filepath, "r")
	lines = f.readlines()
	f.close()
	#print len(lines)
	index = 0
	calcs = []
	for i in range(len(lines)):
		chain=lines[i].split()
		chain.append("XX")
		if chain[1] == "configurations":
			ns = int(chain[0])
			#print ns 
			for j in range(i+1,i+ns+1):
				chain=lines[j].split()
				index = index + 1
				duo = ["Calc"+str(index),int(chain[1])]
				calcs.append(duo)
				#print duo 
			i = i + ns 	
	return calcs

def getEnergy(name):
	os.chdir(name+"/Opt")
	moved = os.path.isfile('./VacMoved')
	if moved:
		#print os.getcwd()
		os.chdir("../..")
		return False,0.0
	else:
		#print os.getcwd()
		if not os.path.isfile("OUTCAR"):
			print("File path {}/Opt/OUTCAR does not exist. Exiting...".format(name))
			sys.exit()
		f = open("OUTCAR", "r")
		lines = f.readlines()
		f.close()
		for line in reversed(lines):
			if "TOTEN" in line:
				chain=line.split()
				energy =  float(chain[4])
				break
		os.chdir("../..")
		return True,energy

def getSpectra(name):
	print "Working in", name
	os.chdir(name+"/DOS")
	if not os.path.isfile("vasprun.xml"):
		print("File path {}/DOS/vasprun.xml  does not exist. Exiting...".format(name))
		sys.exit()
	raw	 = ET.parse('vasprun.xml')
	root = raw.getroot()
	spectraR=[]
	spectraI=[]
	for imag in root.findall("./calculation/dielectricfunction/imag"):
		for pp in imag.findall("array/set/r"):
			kk=pp.text.split()
			spectraI.append([float(kk[0]),(float(kk[1])+float(kk[2])+float(kk[3]))/3.0])
		break
	for imag in root.findall("./calculation/dielectricfunction/real"):
		for pp in imag.findall("array/set/r"):
			kk=pp.text.split()
			spectraR.append([float(kk[0]),(float(kk[1])+float(kk[2])+float(kk[3]))/3.0])
		break
	os.chdir("../..")
	if not spectraR:
  		print("WARNING: The xml file does not contain the dielectric function.")
	if not spectraI:
  		print("WARNING: The xml file does not contain the dielectric function.")
	return spectraR, spectraI

def printspectra(label,spectraR,spectraI):
	c = 29979245800 # cm/s
	planck = 4.135667E-15 # eV/s
	spectraT=[]
	for vecI, vecR in zip(spectraI, spectraR):
		#Refractive index
		refract = math.sqrt(math.sqrt(math.pow(vecR[1], 2)+math.pow(vecI[1], 2))+vecR[1])/math.sqrt(2)
		#Extinction coefficient
		extin =   math.sqrt(math.sqrt(math.pow(vecR[1], 2)+math.pow(vecI[1], 2))-vecR[1])/math.sqrt(2)
		#Absorption coefficient
		absc = 4*math.pi*vecR[0]*extin/(c*planck)
		#Energy loss spectrum
		els = vecI[1]/(math.pow(vecR[1], 2)+math.pow(vecI[1], 2))
		#Reflectivity
		epsilon = complex(vecR[1],vecI[1])
		sqrt_epsi = cmath.sqrt(epsilon)
		reflex = math.pow(abs((sqrt_epsi-1)/(sqrt_epsi+1)),2)
		spectraT.append([vecI[0],vecR[1],vecI[1],refract,extin,absc,els,reflex])

	#Checking f-SumRule
	x = []
	y = []
	for i in spectraT:
		x.append(i[0])
		if i[0] == 0:
                	y.append(0)
		else:
			y.append(i[6]/i[0])
	xx = np.array(x)
	yy = np.array(y)
	integral = trapz(yy, xx)
	print math.pi/2, integral


	filename = label
	f= open(filename,"w+")
	f.write("# Properties extracted from {} using OpticsNML\n".format(label))
	if abs(100*((math.pi/2)-integral)/(math.pi/2))>10:
        	f.write('#WARNING. Check your calculation settings. Sum-rule does not look good: {:>14.4f}\n'.format(integral))
	f.write("# Columns: \n")
	f.write("# 1. Energy (eV) \n")
	f.write("# 2. Dielectric function R \n")
	f.write("# 3. Dielectric function I \n")
	f.write("# 4. Refractive index  \n")
	f.write("# 5. Extinction coefficient  \n")
	f.write("# 6. Absorption coefficient (cm-1)  \n")
	f.write("# 7. Energy loss spectrum  \n")
	f.write("# 8. Reflectivity  \n")
	for i in spectraT:
        	f.write('{:>14.4f}{:>14.4f}{:>14.4f}{:>14.4f}{:>14.4f}{:>18.4E}{:>14.6f}{:>14.6f}{:>18.4E}{:>18.4E}\n'.format(i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7]))
	f.close()

def getDosTotal(name):
	print "Working in", name
	os.chdir(name+"/DOS")
	if not os.path.isfile("vasprun.xml"):
		print("File path {}/DOS/vasprun.xml does not exist. Exiting...".format(name))
		sys.exit()
	raw	 = ET.parse('vasprun.xml')
	root = raw.getroot()
	for ll in root.findall("./calculation/dos/i"):
            fermi = ll.text.split()
            efermi = float(fermi[0])
	dos = []
	for spin in root.findall("./calculation/dos/total/array/set/"):
		sdos = []
		for pp in spin.findall("./r"):
			kk=pp.text.split()
			sdos.append([float(kk[0]),(float(kk[1])),float(kk[2])])
		dos.append(sdos)
	os.chdir("../..")
	for spin in dos:
		for point in spin:
			if point[0] > efermi:
				integral = point[2]
	#Calculating number of total electrons and egap
	for fermi in root.findall("./parameters/separator"):
		for elem in fermi.iterfind('i[@name="NELECT"]'):
			 kk=elem.text.split()
			 dummy = float(kk[0])

	nelect = []
	if len(dos) == 1:
		nelect = [float(dummy)]
	else:
		for spin in dos:
			for point in spin:
				if point[0] > efermi:
					kk = point[2]
					nelect=[kk,dummy-kk]
					break
			break
		for point in dos[1]:
			if point[2]>=nelect[1]:
				fermi=[efermi,point[0]]
				break
	gap=[0.0]*len(dos)
	count = 0
	for spin in dos:
		cband = 0 
		for point in spin:
			if point[0] > fermi[count] and point[2] > nelect[count]:
				cband = point[0]
				#print cband
				break
		gap[count] = cband - fermi[count]
		count = count + 1
	
	return efermi, dos, gap 

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("MasterFile", help="File containing all the OUTSOD files.")
	parser.add_argument("-m", "--mode",  help="Property to be calculated: dos, optics or rdf.")
	parser.add_argument("Temp", help="Temperature. Units: Kelvin", type=float)
	parser.add_argument("-b", "--rdf_bin",  help="Bin width for the rdf. Units: Angstrom ", type=float, default=0.1)
	parser.add_argument("-x", "--rdf_max",  help="Maximun distance to be explored by the rdf. Units: Angstrom", type=float, default=15)
	parser.add_argument("-A", "--rdf_indexA",  help="Element 1 index for rdf", type=int, default=0)
	parser.add_argument("-B", "--rdf_indexB",  help="Element 2 index for rdf", type=int, default=0)
	args = parser.parse_args()

	#Inputs	
	filepath = args.MasterFile # Summary file
	mode = args.mode # Mode
	temp = args.Temp # Temperature
	binw = args.rdf_bin # Bin width
	rdfm = args.rdf_max # RDF length
	rdfA = args.rdf_indexA # RDF index 1
	rdfB = args.rdf_indexB # RDF index 1
	
	#Reading SOD file with degeneracy
	calcs = getList(filepath)

	# Getting energies, ground state and total number of total structures
	ground = 0.0
	ground_label = ""
	totalN = 0
	scalcs=[]
	for calc in calcs:
		valid,energy = getEnergy(calc[0])
		if valid:
			#Saving ground state
			if energy < ground:
				ground = energy
				ground_label = calc[0]
			totalN = totalN + calc[1]
			entry=[calc[0],calc[1],energy]
			scalcs.append(entry)
	print "Ground state is ", ground_label, ", energy: ",ground," eV"
	print "Total number of structures: ",totalN

	#Preparing Energy histogram and calculating DeltaE and partition functions
	efile= open("EnergyHistogram.dat","w+")
	efile.write("# Structure energies (eV) \n")
	pfSum = 0.0
	for calc in scalcs:
		calc.append(calc[2]-ground)
		expo = -1*calc[3]/(temp*KB)
		pf = calc[1]*math.exp(expo)
		pfSum = pfSum + pf
		calc.append(pf)
		for i in range(calc[1]):
			efile.write("%s\n" % calc[3])
	efile.close()
	gnuHist()

	#Sorting by energy	
	scalcs.sort(key=lambda k: (k[3]), reverse=False)
	#Printing list and reducing the list a list of strcutures which represent to P=99.99%
	suma = 0.0
	N = 0
	scalcs2=[]
	for calc in scalcs:
		suma = suma + calc[4]/pfSum 
		scalcs2.append(calc)
		N = N + 1
		if suma >= 0.9999:
			break

	#Calculating optic spectra
	if mode == "optics":
		groundSpectraR = []
		groundSpectraI = []
		#Getting ground state spectra
		groundSpectraR, groundSpectraI=getSpectra(ground_label)
		#Printing ground state spectra
                label = "gsSpectra.dat"
                printspectra(label,groundSpectraR,groundSpectraI)
		#Gettin the averaged spectra
		finalSpectraR = []
		finalSpectraI = []
		#Shaping final Spectra with first calculation
		pointsR, pointsI=getSpectra(scalcs2[0][0])
		for pointR, pointI in zip(pointsR, pointsI):
			pointR[0]= pointR[0]/float(N)
			pointR[1]= scalcs2[0][4]*pointR[1]/pfSum
			finalSpectraR.append(pointR)
			pointI[0]= pointI[0]/float(N)
			pointI[1]= scalcs2[0][4]*pointI[1]/pfSum
			finalSpectraI.append(pointI)
		#Add all the other spectras
		for i in range(len(scalcs2)-1):
			pointsR,pointsI=getSpectra(scalcs2[i+1][0])
			for j in range(len(pointsR)):
				finalSpectraR[j][0] = finalSpectraR[j][0] + pointsR[j][0]/float(N)
				finalSpectraR[j][1] = finalSpectraR[j][1] + scalcs2[i+1][4]*pointsR[j][1]/pfSum
				finalSpectraI[j][0] = finalSpectraI[j][0] + pointsI[j][0]/float(N)
				finalSpectraI[j][1] = finalSpectraI[j][1] + scalcs2[i+1][4]*pointsI[j][1]/pfSum
		#Printing averaged spectra
		filename = "avSpectra"+str(temp)+"K.dat"
                printspectra(filename,finalSpectraR,finalSpectraI)
		gnuSpctra()

	#Calculating optic spectra
	if mode == "dos":
		groundDos = []
		#Gettin the averaged spectra
		finalDos = []
		finalgap = []
		#Shaping final DOS with first calculation
		efermi,groundDos,gap=getDosTotal(scalcs2[0][0])
                print "Participation ratio = ", scalcs2[0][4]/pfSum
		ii = 0
		for spin in groundDos:
			print "Egap spin", ii+1, " = ", gap[ii], " eV"
			spinDOS = []
			finalgap.append(gap[ii]*scalcs2[0][4]/pfSum)
			for point in spin:
				point[0]= point[0]/float(N)
				point[1]= scalcs2[0][4]*point[1]/pfSum
				point[2]= scalcs2[0][4]*point[2]/pfSum
				spinDOS.append(point)
			ii = ii+1
			finalDos.append(spinDOS)
		#Add all the other spectras
		for i in range(len(scalcs2)-1):
			energy,calcDos,gap=getDosTotal(scalcs2[i+1][0])
                	print "Participation ratio = ", scalcs2[i+1][4]/pfSum
		#Updating Fermi energy
			if energy > efermi:
				efermi = energy
			ii = 0
			for j in range(len(calcDos)):
				print "Egap spin", ii+1, " = ", gap[ii], " eV"
				finalgap[ii] = finalgap[ii] + scalcs2[i+1][4]*gap[ii]/pfSum
				for k in range(len(calcDos[0])):
					finalDos[j][k][0] = finalDos[j][k][0] + calcDos[j][k][0]/float(N)
					finalDos[j][k][1] = finalDos[j][k][1] + scalcs2[i+1][4]*calcDos[j][k][1]/pfSum
					finalDos[j][k][2] = finalDos[j][k][2] + scalcs2[i+1][4]*calcDos[j][k][2]/pfSum
				ii = ii+1
		#Printing averaged spectra and applying Fermi energy correction
		spinlabel = 1
		for spin in finalDos:
			filename = "avDOS_spin_"+str(spinlabel)+"_"+str(temp)+"K.dat" 
			f= open(filename,"w+")
			f.write("# Averaged DOS, spin {}\n".format(spinlabel))
			f.write('Fermi energy {:>14.4f}\n'.format(finalgap[spinlabel-1]))
			f.write("# Energy (eV), density of states (states/eV) \n")
			for point in spin:
				f.write('{:>14.6f}{:>14.4f}{:>14.4f}\n'.format(point[0]-efermi,point[1],point[2]))
			f.close()
			spinlabel = spinlabel+1
        #Getting ground state spectra
                print "Printing GS"
		efermi_ground,groundDos,groundgap=getDosTotal(ground_label)
		#Fermi energy correction and printing on file ground state total DOS
		spinlabel = 1
		for spin in groundDos:
			print "Egap spin", spinlabel, " = ", gap[spinlabel-1], " eV"
			filename = "gsDOS_spin_"+str(spinlabel)+".dat" 
			f= open(filename,"w+")
			f.write("# Ground state DOS, spin {}\n".format(spinlabel))
			f.write('Fermi energy {:>14.4f}\n'.format(groundgap[spinlabel-1]))
			f.write("# Energy (eV), density of states (states/eV) \n")
			for point in spin:
				f.write('{:>14.6f}{:>14.4f}{:>14.4f}\n'.format(point[0]-efermi,point[1],point[2]))
			f.close()
			spinlabel = spinlabel + 1
		#Printing standard gnuplot file for visualization
		gnuDOS()

	#Calculating rdf
	if mode == "rdf":
		#Preparing bin sizes and range of the crystal we will explore	
		Nbins = rdfm/binw	
		gsRdfFreq = [0.0]*(int(Nbins)+1)
		avRdfFreq = [0.0]*(int(Nbins)+1)
		#Normalizing intensity bin
		intensity = 1.0/binw
		#Calculate number of neighbours cells we have to explore	
		dim= [0]*3
		param= [0]*3
		latt,elements,struct = readPoscar(scalcs2[0][0])
		param[0],param[1], param[2]= abcCalc(latt)		
		print param 
		for i in range(3):
			dim[i] = int(rdfm/(2*param[i]))+1
		print dim
		#Calculating ground state rdf 
		latt,elements,struct = readPoscar(ground_label)
		distlist = []
		distlist = getRdfList(latt,elements,struct,rdfA,rdfB,dim)
		#Printing average rdf
		filename = "gsRDF_A_"+str(rdfA)+"_B_"+str(rdfB)+".dat" 
		f= open(filename,"w+")
		f.write("# Ground state bond distance list\n")
		f.write("# Distance (A)\n")
		print type(distlist)
		print len(distlist)
		for i in distlist:
			f.write('{:>14.6f}\n'.format(i))
		f.close()
		#Calculating averaged rdf
		base = float(elements[rdfA])
		if rdfA == rdfB:
			prefactor = base/2.0
		else:
			prefactor = base
		for calc in scalcs2:
			freq = [0.0]*(int(Nbins)+1)
			latt,elements,struct = readPoscar(calc[0])
			distlist = getRdfList(latt,elements,struct,rdfA,rdfB,dim)
			print "In calc ",calc[0]
			prob = calc[4]/(pfSum*prefactor)
			print "Probability",calc[1],calc[2],calc[3],calc[4],pfSum,prob 
			freq = getCalcbins(distlist,prob,binw,rdfm)
			avRdfFreq = [sum(x) for x in zip(avRdfFreq, freq)]				
		#Printing average rdf
		filename = "avRDF_A_"+str(rdfA)+"_B_"+str(rdfB)+"_"+str(temp)+"K.dat" 
		f= open(filename,"w+")
		f.write("# Averaged rdf\n")
		f.write("# Distance (A), rho \n")
		xlist = []
		for i in range(len(avRdfFreq)):
			x = (i*binw)+binw
			xlist.append(x)
			xnp = np.array(xlist[:i+1])
			ynp = np.array(avRdfFreq[:i+1])
			integral = trapz(ynp, xnp)
			f.write('{:>14.6f}{:>14.6f}{:>14.6f}\n'.format(x,avRdfFreq[i],integral))
			#f.write('{:>14.6f}{:>14.6f}\n'.format(x,avRdfFreq[i]))
		f.close()
		
			


if __name__ == '__main__':
    main()

""" Porting of MATLAB files to Python 
Astrit Zhushi 2011 a.zhushi@cs.ucl.ac.uk
"""

""" 
   WARRNING: implementation in this file is buggy
   Use iwln_struct.py instead of this 
"""
import numpy as np
import scipy as sp

import struct, subprocess
from decimal import Decimal

from ctypes import *
from pylab import *


# Number of angles differs based on the number of Tx and Rx 
# See Table 7-25i of 802.11n HT standard
num_angels = {"2x1":2, "2x2":2, "1x3":4, "3x2":6, "3x3":6, "4x1":6, "4x2":10, "4x3":12, "4x4":12}

# Contains info such as RSSI per receiver, noise etc
class ChannelInfo :
	def __init__(self, raw_data) :
		# remove NETLINK header stuff
		self.unpacked = raw_data[38:]
		# exract noise a b c, bfee_count each held on an unsigned char
		# 0 - 4
		tmp = struct.unpack("BBBBB", self.unpacked[:5])
		self.noise_a = tmp[0] # 0
		self.noise_b = tmp[1] # 1
		self.noise_c = tmp[2] # 2
		self.bfee_count = tmp[3] + (tmp[4]<<8) # 3

		# extract Nrx, Ntx, rssi_a up to antenna_sel
		# 7 - 19
		tmp = struct.unpack("BBBBBbBBBBBB", self.unpacked[7:19])

		self.Nrx = tmp[0] 		# 7
		self.Ntx = tmp[1] 		# 8
		self.rssi_a = tmp[2]		# 9
		self.rssi_b = tmp[3]		# 10
		self.rssi_c = tmp[4]		# 11
		self.noise = tmp[5]		# 12
		self.agc = tmp[6]			# 13
		self.antenna_sel = tmp[7]					# 14
		self.length = tmp[8] + (tmp[9] << 8)		# 15-16
		self.rate = tmp[10] + (tmp[11] << 8)		# 16-17
		self.raw = self.unpacked[19:]
		# number of subcarriers
		self.nrs = 30
		self.perm = []
		self.perm.append(((self.antenna_sel) & 0x3) + 1)
		self.perm.append(((self.antenna_sel >> 2) & 0x3) + 1)
		self.perm.append(((self.antenna_sel >> 4) & 0x3) + 1)
		#print self.perm
	
	def __str__(self) :
		return "[%d] NOISE(A,B,C)=[%d %d %d] Nrx=%d Ntx=%d RSSI(A,B,C)=[%d %d %d] Noise=%d AGC=%d" % (self.bfee_count, self.noise_a, self.noise_b, self.noise_c, self.Nrx, self.Ntx, self.rssi_a,self.rssi_b, self.rssi_c, self.noise, self.agc)

# Class encapsulating Channel State Information
class CsiUtil :

	def parse(self, ci) :
		index = 0
		remainder = 0
		csi = {}

		for i in range(0, ci.nrs) :

			index = index + 3
			remainder = (index%8)

			Hx = numpy.matrix(numpy.zeros((ci.Nrx, ci.Ntx), complex))

			for r in range(0, ci.Nrx):
				for t in range(0, ci.Ntx) :
					first = (struct.unpack("B", ci.raw[index/8])[0]>>remainder)
					second = (struct.unpack("B", ci.raw[index/8+1])[0]<<(8-remainder))
					#first = (first>>remainder)
					#second = (second<<(8-remainder))
					#print c_double((float)(c_byte((first|second)).value))
					tmp = (c_byte(first|second).value)
					#(struct.unpack("B",payload[index/8])[0] >> remainder) | (struct.unpack("B",payload[index/8+1])[0] << (8-remainder))
					real = (c_double(tmp).value)
					first = (struct.unpack("B",ci.raw[index/8+1])[0] >> remainder) 
					second = (struct.unpack("B",ci.raw[index/8+2])[0] << (8-remainder))
					tmp = (c_byte(first|second).value)
					imag = (c_double(tmp).value)
					index+=16
					Hx[r,t] = complex(real,imag)
			csi[str(i)] = Hx

		return csi

	def add(self, subcarrier, cHx):
		key = str(subcarrier)
		self.Hx[key] = cHx

	# print CSI 
	def printCsi(self, ci) :
		for k in range (0, ci.nrs) :
			print "Group %d - #TX=%d #RX=%d bfee_count=%d RSSI A=%d B=%d C=%d NOISE=%d AGC=%d RATE=%d" % (k, ci.Ntx, ci.Nrx, ci.bfee_count, ci.rssi_a, ci.rssi_b, ci.rssi_c, ci.noise, ci.agc, ci.rate)
			print self.parse(ci)#Hx[str(k)]);

	def getPlotData(self, ci) :
		data = []
		#Hx = self.getScaledCsi(self.parse(ci), ci)
	
		Hx = self.newGetScaledCsi(self.parse(ci),ci);

		#Hx = self.parse(ci)
		for i in range(0,ci.nrs) :
		 	key = str(i)
			tmp = abs(Hx[key])
			#tmp = Hx[key]
			#print tmp
			for r in range(0, ci.Nrx) :
				for t in range(0, ci.Ntx ) :
					tmp[t,r] = tmp[t,r] * tmp[t,r]
			d = (10*numpy.log10(tmp))
			#print d.item(0);
			data.append(d.item(0))
		return data

	def getPinv(self, csi, nrs) :
		Hx = csi #getScaledCsi()
		for k in range (0, nrs) :
			key = str(k)
			Hx[key] = numpy.linalg.pinv(Hx[key])
		return Hx
	
	def getScaledCsi(self, Hx, ci) :
		csi_sq = {}
		csi_mag = {}
		scale = {}
		ret = {}
		rssi_mag = 0
		if not ci.rssi_a == 0 :
			rssi_mag = rssi_mag + self.dbinv(ci.rssi_a)
		if not ci.rssi_b == 0 :
			rssi_mag = rssi_mag + self.dbinv(ci.rssi_b)
		if not ci.rssi_c == 0 :
			rssi_mag = rssi_mag + self.dbinv(ci.rssi_c)
		
		noise = 0
		
		if ci.noise == -127 :
			noise = -92
		else :
			noise = ci.noise
		
		for i in range(0,ci.nrs) :
			key = str(i)
			csi_sq[key] = numpy.multiply(Hx[key], numpy.matrix.conj(Hx[key]))
			csi_mag[key] = numpy.matrix.sum(csi_sq[key])
			scale[key] = rssi_mag / (csi_mag[key]/ci.nrs) / self.dbinv(noise+44+ci.agc)
			ret[key] = numpy.multiply(Hx[key], numpy.sqrt(scale[key]))
			if ci.Ntx == 2 :
				ret[key] = numpy.multiply(ret[key], sqrt(ci.Ntx))
			elif ci.Ntx == 3 :
				ret[key] = numpy.multiply(ret[key],sqrt(self.dbinv(4.5)))
		return ret

	#### verify me ####	
	def newGetScaledCsi(self, Hx, ci) :
		rssi = []
		noise = []
 
		csi_sum = 0

		if ci.rssi_a > 0 :
			rssi.append(ci.rssi_a)
		if ci.rssi_b > 0 : 
			rssi.append(ci.rssi_b)
		if ci.rssi_c > 0 :
			rssi.append(ci.rssi_c)

		rssi = [ r - 44 - ci.agc for r in rssi ]

		noise.append(ci.noise_a)
		noise.append(ci.noise_b)
		noise.append(ci.noise_c)
		
		noise = [ noise[i] for i in range(0, ci.Nrx) ]

		ref_noise = max(noise)

		noise_diff = [ ref_noise - n for n in noise ]
		noise_diff_abs = [ numpy.power(10, n/10) for n in noise_diff ]
		ref_rssi = [ r - ref_noise for r in rssi ]

		rssi_sum = numpy.sum( [ pow(10, (r/10)) for r in ref_rssi ] )
		#rssi_sum = numpy.sum(rssi_sum)

		for i in range(0, ci.nrs) :
			key = str(i)
			tmpAbs = numpy.abs(Hx[key])
			csi_sum = csi_sum + numpy.sqrt(tmpAbs)

		common_scale = numpy.sqrt(rssi_sum/csi_sum * ci.nrs)
		scale_per_rx = [ common_scale * numpy.sqrt(n) for n in noise_diff_abs ]

		ret = {}

		for j in range(0, ci.nrs) :
			key = str(j)
			tmpHx = Hx[key]
			tmpRet =  numpy.matrix(numpy.zeros((ci.Ntx, ci.Nrx), complex))
			for r in range(0, len(scale_per_rx)) :
				for t in range(0,ci.Ntx ) :
					tmpVal = tmpHx.item((r,t)) * scale_per_rx[r]
					tmpRet.itemset((t,r), tmpVal.item(0))
			ret[key] = tmpRet

		return ret;
		
	
	def dbinv(self, x) :
		p = numpy.power(10,(x/float(10)))
		return p


   def compress_v(self, csi, psi_bits=3) :
      
	def cal_quantized_steering_matrix(self, csi, phi_bits=5, psi_bits=3) :
		global num_angels
		Nrx = 3#ci.Nrx
		Ntx = 3#ci.Ntx
		key = str(Ntx)+"x"+str(Nrx)
		# number of subcarriers
		nsc = 30#ci.nrs
		# calculate pseudo-inverse
		Hx = csi#self.newGetScaledCsi(csi,ci)#csi#self.getPinv(csi,nsc)
		#for i in range(0,nsc) :
		#	s = str(i)
	#		Hx[s] = Hx[s].getT();

		#Hx = self.getScaledCsi()
		
		ks = numpy.matrix(numpy.zeros((nsc, num_angels[key])))
		bitpattern = numpy.matrix(numpy.zeros((nsc,1)))
		pattern_len = (phi_bits+psi_bits)*3
		bit_spec=[]

		if key == "2x1" or key == "2x2" :
			bit_spec = [psi_bits, phi_bits]
		elif key == "3x1" or key == "1x3":
			bit_spec = [psi_bits, psi_bits, phi_bits, phi_bits]
		elif key == "3x2" or key == "3x3" :
			bit_spec = [psi_bits, phi_bits, psi_bits, psi_bits, phi_bits, phi_bits]
		else :
			print "Only 2x1,2x2,3x1,3x2,3x3 configuration is supported. Give %s" % key

		#print bit_spec
		vcomps = []
		#vcomps.append("../../tools/send_vcomp")
		for i in range(0, nsc) :
			k = str(i)
			(angles, quant_angles, quant_bits) = self.calc_compressed_quantized_angles(Hx[k], phi_bits, psi_bits, Ntx, Nrx)
			quant_bits = [ int(x) for x in quant_bits ]
			#print quant_bits
			#quant_bits.reverse()
			#print quant_bits
			vcomps.append(self.reverse_bits(self.concat_bits(quant_bits, bit_spec), pattern_len))

		#args = " ".join(["%d" % d for d in vcomps])
		#subprocess.call(vcomps)
		return vcomps

	def calc_compressed_quantized_angles(self, m, phi_bits, psi_bits, Ntx, Nrx) :
		U,s,V = numpy.linalg.svd(m)
	
		V = V.T

		theta = numpy.angle(V)
		#theta = numpy.angle(m.T)
		Dtilde = numpy.ma.identity(3,complex)
		
		quant_bits = []
		### see 802.11n Table 7-25i Section 7.3.1.29
		if Ntx == 3 and Nrx == 3:
			quant_bits = [0,0,0,0,0,0]
		elif Ntx == 3 and Nrx == 1 :
			quant_bits = [0,0,0,0]
		elif Ntx == 2 and Nrx == 1 :
			quant_bits = [0,0]

		phi_11 = -(theta[Ntx-1,0]-theta[0,0])
		
		if phi_11 < 0 :
			phi_11 = phi_11 + 2 * math.pi
	
		(phi_11_quant, quant_bits[0]) = self.quantize_phi(phi_11, phi_bits)

		phi_21 = -(theta[Ntx-1,0]-theta[1,0])

		if phi_21 < 0 :
			phi_21 = phi_21 + 2 * math.pi
		
		(phi_21_quant, quant_bits[1]) = self.quantize_phi(phi_21, phi_bits)
		
		d1 = numpy.diag(array([numpy.exp(+1j*phi_11), numpy.exp(+1j*phi_21), 1]), k=0);
		new_v = numpy.matrix.conj(d1) * V

		for i in range(0,len(new_v)):
			new_v[i,0]  = new_v[i,0] * abs(V[2,0])/V[2,0]

		Dtilde[0,0] = V[Ntx-1,0]/abs(V[Ntx-1,0])
		denom = new_v[0, 0].real
		
		if denom != 0 :
			psi_21 = numpy.arctan(new_v[1, 0].real/denom);
		else :
			psi_21 = 0;

		(psi_21_quant, quant_bits[2]) = self.quantize_psi(psi_21, psi_bits)

		G21 = numpy.matrix([[numpy.cos(psi_21),numpy.sin(psi_21),0],[-numpy.sin(psi_21),numpy.cos(psi_21),0],[0,0,1]])
		v21 = G21 * new_v;
		psi_31 = numpy.arctan(v21[2, 0].real/v21[0, 0].real);
		
		(psi_31_quant, quant_bits[3]) = self.quantize_psi(psi_31, psi_bits)
		
		G31 = numpy.matrix([[numpy.cos(psi_31),0, numpy.sin(psi_31)],[0,1,0],[-sin(psi_31),0,cos(psi_31)]]);
		v31 = G31 * v21;
			
		# phi_22
		theta_32 = numpy.angle(v31[2,1]);
		theta_22 = angle(v31[1, 1]);

		phi_22 = -(theta_32 - theta_22);
		if phi_22 < 0 :
			phi_22 = phi_22 + 2 * math.pi

		(phi_22_quant, quant_bits[4]) = self.quantize_phi(phi_22, phi_bits)

		d2 = numpy.diag(array([1, numpy.exp(+1j*phi_22),1]));

		v22 = numpy.matrix.conj(d2) * v31;

		for i in range(0, len(v22)) :
			v22[i, 1] = v22[i, 1] * abs(v31[2, 1])/v31[2,1]

		Dtilde[1, 1] = v31[2,1]/abs(v31[2,1])

		psi_32 = numpy.arctan(v22[2,1].real/v22[1, 1].real)

		(psi_32_quant, quant_bits[5]) = self.quantize_psi(psi_32, psi_bits)
			
		G32 = numpy.matrix([[1,0,0],[0,numpy.cos(psi_32), numpy.sin(psi_32)],[0,-numpy.sin(psi_32),numpy.cos(psi_32)]])
		v32 = G32 * v22;
		Dtilde[2, 2] = v32[2, 2] / abs(v32[2, 2])

		angles = [psi_21, psi_31, phi_11, phi_21]
		#angles = [ phi_11, psi_21]
		angles = [phi_11,phi_21,psi_21,psi_31,phi_22,psi_32];
		quant_angles = [phi_11_quant,phi_21_quant,psi_21_quant,psi_31_quant,phi_22_quant,psi_32_quant];
		#quant_angles = [ phi_11_quant,psi_21_quant];
		return (angles,quant_angles, quant_bits)

	# quantize phi according to 7.3.1.29, assuming phi within [0, 2pi] 
	def quantize_phi(self, phi, phi_bits) :
		kmax = 2**phi_bits - 1;
		phi_min = pi/(2**phi_bits);
		phi_max = kmax*math.pi/2**(phi_bits - 1) + phi_min;

		if phi <= phi_min:
			return (phi_min, 0)

		if phi >= phi_max:
			return (phi_max, kmax)

		unit = pi / (2**(phi_bits - 1))

		num_units = phi / unit

		k = floor(num_units)
 
		return (unit * (k + 0.5), k);
			

	# quantize psi according to 7.3.1.29, assuming phi within [0, 2pi] 
	def quantize_psi(self, psi, psi_bits) :
		kmax = numpy.power(2,psi_bits)-1;
		psi_min = math.pi/(numpy.power(2,(psi_bits+2)));
		psi_max = kmax*math.pi/numpy.power(2,(psi_bits + 1)) + psi_min;
		
		if psi <= psi_min:
			return (psi_min, 0)

		if psi >= psi_max:
			return (psi_max, kmax)

		unit = pi / (numpy.power(2,(psi_bits + 1)))

		num_units = psi / unit
		k = numpy.floor(num_units)
		unit = unit * (k+0.5)
		ret = (unit, k)
		return ret


	def concat_bits(self, numbers, sizes) :
		#
		# Given a list of numbers and how many bits they each take, 
		# concatenate them into a long bit string.
		#
		length = len(numbers)

		if length != len(sizes) :
			print 'Lengths of input arguments don''t match!'
			return 0
		
		res = numbers[0];
		for idx in range(2,length) :
			res = (res<<sizes[idx])
			res = res + numbers[idx]
		return res

	def reverse_bits(self, orig, num_bits) :
		#
		# reverses the bits in orig from msb->lsb to lsb->msb, 
		# assuming orig has num_bits.
		#
		min_bits = numpy.floor(numpy.log2(orig)) + 1;
		max_bits = 32;
		res = 0;

		if num_bits < min_bits :
			print 'Need at least %d bits' % min_bits
			return
		elif num_bits > max_bits :
			print 'Only supports numbers up to %d bits!' % max_bits
			return

		temp1 = ( orig & int('0xaaaaaaaa',0)) >> 1
		temp2 = ( orig & int('0x55555555',0)) << 1
		orig =  (temp1 | temp2)

		temp1 = (orig & int('0xcccccccc',0)) >> 2
		temp2 = (orig & int('0x33333333',0)) << 2
		orig =  (temp1 | temp2)

		temp1 = (orig & int('0xf0f0f0f0',0)) >> 4
		temp2 = (orig & int('0x0f0f0f0f',0)) << 4
		orig = (temp1 | temp2)

		temp1 = (orig & int('0xff00ff00',0)) >> 8
		temp2 = (orig & int('0x00ff00ff',0)) << 8
		orig = (temp1 | temp2)

		temp1 = orig >> 16
		temp2 = orig << 16 
		
		orig = (temp1 | temp2)

		bitmask = 2**num_bits - 1;
		shift = int(numpy.ceil(num_bits)) - max_bits

		if ( shift > 0 ) :
			temp = orig << shift
		else :
			temp = orig >> (-1*shift)

		r = (temp & bitmask)
		return r

	def printMatlabReady(self, m) :
		print "[",
		for r in range(0,self.ci.Nrx) :
			for c in range(0, self.ci.Ntx) :
				print m[r,c],
			print ';',
		print "]",
			


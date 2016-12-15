#!/usr/bin/python

import socket,os
import platform

""" NETLINK related stuff
Astrit Zhushi 2011, a.zhushi@cs.ucl.ac.uk
"""

NETLINK_CONNECTOR=11
NETLINK_ADD_MEMBERSHIP=1

def get_cn_idx_iwlagn():
	uname = platform.uname()[2]
	infile = open("/usr/src/linux-headers-%s/include/linux/connector.h"
			%(uname), "r")
	flag = False
	for line in infile:
		if line.find("CN_IDX_IWLAGN") == -1:
			continue
		line = line.strip().split()
		CN_IDX_IWLAGN = eval(line[2])
		flag = True
		break
	infile.close()

	if flag:
		return CN_IDX_IWLAGN

	raise IOError("CN_IDX_IWLAGN not found in connector.h")

def get_iwlnl_socket() :
	CN_IDX_IWLAGN = get_cn_idx_iwlagn()
	s = socket.socket(socket.AF_NETLINK, socket.SOCK_DGRAM, NETLINK_CONNECTOR)
	pid = os.getpid()
	s.bind((pid,CN_IDX_IWLAGN))
	s.setsockopt(270, NETLINK_ADD_MEMBERSHIP, CN_IDX_IWLAGN)
	return s



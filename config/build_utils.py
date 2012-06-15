#! /usr/bin/env python
################################################################################
# Copyright (c) 2012 Los Alamos National Security, LLC
# All rights reserved.
#
# $Revision: 43 $
# $Date: 2012-03-27 20:49:25 -0600 (Tue, 27 Mar 2012) $
# $Author: bergen $
################################################################################

import os, string, sys, shutil

class Test:
	'''
	Test executable and the source files that build it.
	'''
	def __init__(self, exe, srcs, defines, procs):
		self.exe = exe
		self.srcs = string.split(srcs)
		self.defines = defines
		self.procs = procs

class Test_Builder:
	'''
	Class for building tests
	'''
	def __init__(self):
		# regular tests
		self._have_test = False
		self._tests = []
		self._sources = []
		# input files and dirs
		self._test_input_files = []
		self._test_input_dirs = []

	def test(self, exe, srcs, defines, procs=[1]):
		self._tests.append(Test(exe, srcs, defines, procs))
		self._have_test = True

	def test_input_files(self,files):
		self._test_input_files.extend(string.split(files))

def get_bld(file, bld):
	'''
	Parses file for bld.* statements.
	'''
	execfile(file)
	return bld

if __name__ == '__main__':
	_usage ='''
	Usage: config/build_utils.py CONFIGDIR
	'''

	if len(sys.argv) != 2 :
		print _usage
		sys.exit(1)

	configdir = sys.argv[1]

	topdir = os.getcwd()

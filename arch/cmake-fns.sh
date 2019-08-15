#!/usr/bin/env bash
#
# Copyright (c) 2019 Carnegie Mellon University,
# Copyright (c) 2019 Triad National Security, LLC, as operator of
#     Los Alamos National Laboratory.
#
# All rights reserved.
#
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file. See the AUTHORS file for names of contributors.
#

#
# cmake-fns.sh  cmake helper functions
# 18-Jul-2019  chuck@ece.cmu.edu
#

#
# create_cmake_cache: create a new cmake cache
#
create_cmake_cache() {
    init_cache=./vpic-init-cache.cmake
    rm -f $init_cache
    touch $init_cache
    if [ ! -w $init_cache ]; then
        echo "ERROR: failed to create $init_cache init cache"
        exit 1
    fi
}

#
# set_cache_bool: add a bool to the current cache
#
set_cache_bool() {
  echo "set($1 \"$2\" CACHE BOOL \"Initial cache\" FORCE)" >> $init_cache
}

#
# set_cache_string: add a string to the current cache
#
set_cache_string() {
  echo "set($1 \"$2\" CACHE STRING \"Initial cache\" FORCE)" >> $init_cache
}

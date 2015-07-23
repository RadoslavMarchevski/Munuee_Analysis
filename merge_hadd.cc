#!/bin/bash

tag=$1
version=$2

hadd $1_merged_$2.root $1*/*.root

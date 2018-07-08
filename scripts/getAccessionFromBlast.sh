#!/bin/bash

grep -i access $1 | grep -o '>.*' | grep -o '[^>]*' | grep -o '^[^<]*'

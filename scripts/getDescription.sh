#!/bin/bash

grep -o '[^>]*' $1 | cut -d ' ' -f 2- > $2

#!/bin/bash

sed -r 's/^(>[A-Za-z_0-9]+)\s.+$/\1/g' Asax_cdhit.fa > Asax_truncnames.fa

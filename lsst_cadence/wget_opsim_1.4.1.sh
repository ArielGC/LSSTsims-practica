#!/bin/bash

wget -np -nc -nd https://epyc.astro.washington.edu/~lynnej/opsim_downloads/fbs_1.4/ -r -k -e robots=off -p -B https://epyc.astro.washington.edu/~lynnej/opsim_downloads/fbs_1.4/ -A '*10yrs.db' -P ./FBS_1.4.1/
#!/bin/bash

$input_fa=$1
$db_fa=$2
$output=$3

igblastp -query $input_fa -germline_db_V $db_fa -out $output -outfmt 7 -num_threads 24
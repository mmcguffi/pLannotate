#!/bin/bash
while IFS= read -r line; do
  name=`echo -n $line | sed "s|https://media.addgene.org/snapgene-media/v1.6.2-0-g4b4ed87/sequences||"`
  name=`echo -n $name | sed "s|/|_|g"`
  #echo wget "$line" -P ./addgene_gbks/
  curl "$line" -so ./addgene_gbks/$name
done < "$1"
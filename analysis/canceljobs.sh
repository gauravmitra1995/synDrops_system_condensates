#!/bin/bash

for id in `seq 40097819 1 40097833`;do
    scancel $id
done


#!/bin/bash

cd /yljobs

find . -mtime +14 -type f -delete
find . -mtime +14 -type d -exec rm -Rf {} +

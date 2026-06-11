#!/bin/bash

# Ensure we have the latest info from the server
git fetch --all --quiet

# Define Header
printf "%-25s | %-20s | %-15s | %-20s | %-20s\n" "Branch" "Creator" "Date Created" "Last Committer" "Relative Age"
printf "%0.s-" {1..105}
printf "\n"

# Use for-each-ref to get remote branches sorted by creation date (oldest first)
# 'creatordate' refers to the date the branch was first pushed/created
for ref in $(git for-each-ref --sort=creatordate --format='%(refname:short)' refs/remotes/origin | grep -v 'HEAD'); do
    
    # Remove 'origin/' prefix for display
    branch_name=$(echo $ref | sed 's/origin\///')
    
    # 1. Get the Creator (First commit on this branch not in master)
    creator=$(git log master..$ref --reverse --format="%an" | head -n 1)
    created_date=$(git log master..$ref --reverse --format="%ad" --date=short | head -n 1)
    
    # 2. Get the Latest Info
    last_committer=$(git log -1 $ref --format="%an")
    relative_date=$(git log -1 $ref --format="%ar")

    # Handle cases where the branch is already merged or identical to master
    if [ -z "$creator" ]; then
        creator="Merged/Same"
        created_date="N/A"
    fi

    # Print aligned row
    printf "%-25s | %-20s | %-15s | %-20s | %-20s\n" "$branch_name" "$creator" "$created_date" "$last_committer" "$relative_date"
done
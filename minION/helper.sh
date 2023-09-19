#!/bin/bash
# helper.sh

# This is a helper function to join two paths
join_paths() {
    local path1="$1"
    local path2="$2"

    # Join paths
    local joined_path="${path1}/${path2}"

    # Normalize path (removing double slashes if they exist)
    joined_path=$(echo "$joined_path" | sed 's#//+#/#g')

    echo "$joined_path"
}
cd "$(dirname "$0")" || exit 1
rsync -avz --progress ./sync GPU-PC:~/

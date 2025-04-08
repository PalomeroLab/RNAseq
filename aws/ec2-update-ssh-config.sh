#!/bin/sh

# SSH Config File
SSH_CONFIG="$HOME/.ssh/config"

# Function to update the Hostname in the SSH config file
update_ssh_config() {
	# Check if the host exists in the SSH config file
	if grep -q "Host $host" "$SSH_CONFIG"; then
		# If the host exists, replace the Hostname
		sed -i "s/^  Hostname .*/  Hostname $new_hostname/" "$SSH_CONFIG"
		echo "Updated Hostname for $host to $new_hostname."
	else
		# If the host doesn't exist, print a message
		echo "Host $host not found in the SSH config file."
	fi
}

# Get the directory of the current script
this_dir=$(dirname "$(readlink -f "$0")")

# Change to the directory where the script is located
cd "$this_dir" || {
	printf "Failed to change directory to %s\n" "$this_dir" >&2
	exit 1
}

# Select the EC2 instance and get the relevant info (host and hostname)
selected_instance=$(../aws/ec2-list-instances.sh | fzf | awk '{if (NF >= 3) print $1, $3; else print $1, "No Hostname"}')

# Check if a valid selection was made
if [ "$selected_instance" = "" ]; then
	printf "No instance selected. Exiting.\n" >&2
	exit 1
fi

# Capture the host and hostname from the selected line
host=$(echo "$selected_instance" | awk '{print $1}')
hostname=$(echo "$selected_instance" | awk '{print $2}')

# Check if hostname is valid
if [ "$hostname" != "No Hostname" ]; then
	# Update the SSH config with the new hostname
	echo "Selected instance: Host: $host, Hostname: $hostname"
	update_ssh_config
else
	echo "No Hostname found for instance $host, skipping SSH config update."
fi

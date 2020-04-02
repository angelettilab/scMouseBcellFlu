#! /bin/bash

# To build remotely, you need to first connect to a sylabs remote server using a token.
# Run the following command:
#
#singularity remote login
#
# Then paste the API key (in the "sylabs-token" file) into the terminal when prompted.
# You should see the message "API Key Verified!".

# build image
singularity build --remote sauron.sif Singularity_remote.def

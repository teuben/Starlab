{
    if ($2 == option && $1 != "unsetenv") {
	if (action > 0) 
	    printf "setenv\t";
	else
	    printf "#setenv\t";

#	/* Note excruciating formatting! */

	printf "%s\t", $2;
	for (i = 3; i <= NF; i++) {
	    if (i == NF)
		print $NF;
	    else if ($i == "#") {
		if (i == 4) printf "\t";
		printf "\t# ";
	    } else
		printf "%s ", $i;
	}

    } else
	print $0;
}

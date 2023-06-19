for vib_mode in {58..72..1}; do
	for vib_amp in $(seq -1.5 0.5 1.5); do
		cp mesh.conf mesh_amp.conf;
		echo $vib_freq $amp
		sed -i "s/vib_mode/$vib_mode/g" mesh_amp.conf
		sed -i "s/vib_amp/$vib_amp/g" mesh_amp.conf
		tail -n 1 mesh_amp.conf
		phonopy mesh_amp.conf;
		cp MPOSCAR-orig MPOSCARorig$vib_mode$vib_amp;
		cp MPOSCAR MPOSCAR$vib_mode$vib_amp;
		mv MPOSCARorig$vib_mode$vib_amp /home/armin/Downloads/Exjobb/Basic_Tutorial/MPOSCAR;
		mv MPOSCAR$vib_mode$vib_amp /home/armin/Downloads/Exjobb/Basic_Tutorial/MPOSCAR;
	done
done

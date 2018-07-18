import sys;

# open prosite.dat file and read-in all lines
f_name = sys.argv[1];
f = open(f_name,'r');
lines= [];
line = f.readline();
while line:
	lines.append(line);
	line = f.readline();
f.close();

#define classes
classes = [ "cTP", "cyt","ER","SR","gol","lys","mTP","nuc","per","mem","vac"];

# read-in clusterfile and save all clusters
f_name2 = sys.argv[2];
f = open(f_name2,'r');
line = f.readline();
c = -1;
while line:
	if line[0:11] == "clust_class":
		c += 1;
		#extract PS ids
		line_s = line.split("[");
		line_s = line_s[1].split("]");
		line_s = line_s[0].split(",");
		for i in range(len(line_s)):
			line_s[i] = line_s[i].strip();
		print(line_s);
		#write prosite id file
		f_id_name = "prosite_ids.txt_";
		if c < len(classes):
			f_id_name += "clust_"+classes[c];
		else:
			f_id_name += "anti_clust_"+classes[c-len(classes)];
		f_id = open(f_id_name,'w');
		for elem in line_s:
			f_id.write(elem+"\n");
		f_id.close();
		# set write boolean to True
		write = True;
		#create new dat file
		f_out_name = f_name+"_";
		if c < len(classes):
			f_out_name += "clust_"+classes[c];
		else:
			f_out_name += "anti_clust_"+classes[c-len(classes)];
		f_out = open(f_out_name,'w');
		# iterate through prosite.dat file and check whether entry is written or not
		for i in range(len(lines)):
			if i+1 >= len(lines):
				write = True;
			elif lines[i][0:2] == "//":
				p_id = lines[i+2][5:12];
				if p_id in line_s:
					write = True;
				else:
					write = False;
			if write == True:
				f_out.write(lines[i]);
		f_out.close();
	line = f.readline();
f.close();

/* 
* Copyright (C) <2026>  <The Ohio State University>       
* 
* This program is free software: you can redistribute it and/or modify                              
* it under the terms of the GNU General Public License as published by 
* the Free Software Foundation, either version 3 of the License, or    
* (at your option) any later version.                                                                                       
* This program is distributed in the hope that it will be useful, 
* but WITHOUT ANY WARRANTY; without even the implied warranty of           
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       
* GNU General Public License for more details.                                                                             
* 
* You should have received a copy of the GNU General Public License 
* along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/ 

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <math.h>

extern "C" {
#include <ViennaRNA/model.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/utils/strings.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/params/io.h>
}

using namespace std;

//compile with
//g++ RNAfec.c -I ../src -L ../src/ViennaRNA -o RNAfec.x -lm -lgomp -lRNA -fopenmp

//void usage(void) {
//  cout << "usage: ./RNAfec.x <Sequence> <Protein name> <conc>\n";
//}

struct options {
  double          min_force;
  double          max_force;
  double          delta_force;
  vrna_md_t       md;
  char            *output_file;
};

struct record_data {
  char            *sequence;
  char            *protein_name;
  char            *input_filename;
  double          conc;
  int             num_motifs;
  string          *ligandMotifs;
  double          *ligandinvkds;
  int             num_fs;
  double          *energies;
  double          *forces;
  struct options  *options;
};

void
init_default_options(struct options *opt)
{
  opt->min_force      = 0.0;
  opt->max_force      = 25.0;
  opt->delta_force    = 0.02;
  set_model_details(&(opt->md));
  opt->output_file        = NULL;
}

static void 
compute_Energies (struct record_data *record);

static void 
compute_exten_cd (struct record_data *record);

int
main(int argc, char * argv[])
{ 
  struct options *opt = (struct options *)vrna_alloc(sizeof(struct options));
  unsigned int  count;
  char line[256];                                                                                                                                                      
  double invKd;

  init_default_options(opt);
	
  struct record_data *record = (struct record_data *)vrna_alloc(sizeof(struct record_data));
  
  record->sequence        = argv[1];
  printf(" %s \n", record->sequence);
  fflush(stdout);
  record->options         = opt;
  record->conc            = 0;
  record->num_fs		  = (int) (((opt->max_force - opt->min_force) / opt->delta_force) + 1);
  record->energies        = (double *) vrna_alloc(sizeof(FLT_OR_DBL)*record->num_fs);
  record->forces          = (double *) vrna_alloc(sizeof(FLT_OR_DBL)*record->num_fs);

  if (argc > 2){
	record->protein_name    = argv[2];
	record->conc            = atof(argv[3]);
	if (record->conc > 0){
		char completefilename[256];
		const char path[]             = "./protein_files/";
		const char filetype[]         = ".rnacompete";
		strcpy(completefilename, path);
		strcat(completefilename, record->protein_name);
		strcat(completefilename, filetype);
		record->input_filename  = completefilename;
		printf(" %s\n", record->input_filename);
		fflush(stdout);
		printf(" %f\n", record->conc);
		fflush(stdout);
		record->num_motifs      = 0;
		record->ligandMotifs    = (string *) vrna_alloc(550502400);
		record->ligandinvkds    = (double *) vrna_alloc(550502400);
		
		FILE* lengthlist = fopen("./protein_files/RBDls.txt", "r");
		
		while (fgets(line, sizeof(line), lengthlist)) {
			/* note that fgets don't strip the terminating \n, checking its
		  presence would allow to handle lines longer that sizeof(line) */
		  char* token;
		  char* rest = line;
		  char* protein;
		  int argcount = 0;
	 
		  while ((token = strtok_r(rest, " ", &rest))){
			fflush(stdout);
			if(argcount==0){
				protein = token;
			}
			else if(argcount==1){
			  int comp = strcmp(record->protein_name,protein);
			  if (comp==0){
				opt->md.protein_l = atof(token);
				printf(" %f\n", opt->md.protein_l);
				fflush(stdout);
			  }
			}
			argcount++;
		  }
		}

		FILE* motiflist = fopen(record->input_filename, "r");
		count = 0;
		
		printf("%s\n", record->input_filename);
		fflush(stdout);
	  
		while (fgets(line, sizeof(line), motiflist)) {
		  /* note that fgets don't strip the terminating \n, checking its
		  presence would allow to handle lines longer that sizeof(line) */
		  char* token;
		  char* rest = line;
		  int argcount = 0;
	 
		  while ((token = strtok_r(rest, " ", &rest))){
			if(argcount==0){
				record->ligandMotifs[count] = token;
			}
			else if(argcount==1){
			  char *eptr;
			  invKd = strtod(token, &eptr);
			  record->ligandinvkds[count] = invKd;
			}
			argcount++;
		  }
		
		  count++;
		}
	  
		record->num_motifs = count;
		printf("%d\n", record->num_motifs);
		fflush(stdout);
	}
  }
  
  /* convert sequence to uppercase letters only */
  vrna_seq_toupper(record->sequence);
  
  printf("Calculating Energies... \n");
  fflush(stdout);
  /* Calculate Energies */
  compute_Energies(record);
  
  printf("Calculating Extensions... \n");
  fflush(stdout);
  /* Calculate data for FEC */
  compute_exten_cd(record);
  
  /*clean up*/
  free(opt);
  free(record);

  printf("All done... \n");
  fflush(stdout);

  return EXIT_SUCCESS;

}

static void
compute_Energies(struct record_data *record){
	
	unsigned int          length, f_index, max_f_index;
    struct options        *opt;
    char                  *rec_sequence; 
	string                *motifs;//, mfe_structure;
    double                ens_energy, min_f, delta_f, conc, *invkds;//, min_en;
    vrna_fold_compound_t  *vc;
	
	opt              = record->options;
	min_f            = opt->min_force;
    delta_f          = opt->delta_force;
    max_f_index      = record->num_fs;
	rec_sequence     = strdup(record->sequence);
	motifs 		     = record->ligandMotifs;
	invkds	    	 = record->ligandinvkds;
	conc 			 = record-> conc;

	for(f_index=0; f_index < max_f_index; f_index++) {
		/* calculate force */
		double f = min_f + delta_f * f_index;

		/* update force stored in model details */
		opt->md.force = f;
		/*opt->md.temperature = 24.85;*/
		record->forces[f_index] = f;


		/* create fold compound */	
		vc = vrna_fold_compound(rec_sequence, &(opt->md), VRNA_OPTION_DEFAULT);

		length = vc->length;

		//mfe_structure = (char *)vrna_alloc(sizeof(char) * (length + 1));

		/*
		########################################################
		# begin actual computations
		########################################################
		*/

		//min_en = (double)vrna_mfe(vc, mfe_structure);

		//if (length > 2000)
		//  vrna_mx_mfe_free(vc);

		char *pf_struc = (char *)vrna_alloc(sizeof(char) * (length + 1));

		//vrna_exp_params_rescale(vc, &min_en);


		if (length > 2000){
			vrna_message_info(stderr, "scaling factor %f", vc->exp_params->pf_scale);
		}
		
		if (conc > 0){
			/*add ligand*/
			for(int s=0;s < record->num_motifs;s++) {
				const char *char_motif = motifs[s].c_str();
				double energy = 8.3144598*310.15*log((1.0/invkds[s])/conc)/4184.0;
				vrna_ud_add_motif(vc, char_motif, energy, NULL, 'A');
			}
		}

		/* calculate and store ensemble energy in pN nm w/o ligand */
		ens_energy = (double) vrna_pf(vc, pf_struc)*6.9513533;
		record->energies[f_index] = ens_energy;

		vrna_ud_remove(vc);

		free(pf_struc);

		/* clean up */
		vrna_fold_compound_free(vc);


	}
	
	free(rec_sequence);
}

static void 
compute_exten_cd (struct record_data *record)
{
	unsigned int          num_fs, index, h_ind;
	struct options        *opt;
	double                exten, force, min_f, max_f, delta_f, h;
	double 				 *forces, *energies;

	opt      = record->options;
	forces   = record->forces;
	energies = record -> energies;
	min_f    = opt->min_force;
	max_f	   = opt->max_force;
	delta_f  = opt->delta_force;
	num_fs   = record->num_fs;
	h        = delta_f * 2.;
	
    char filename[16384];
	char concname[318];
	snprintf(concname, sizeof concname, "%f", record->conc);
	char end[] = ".txt";
	char conn[] = "_";
	strcpy(filename, "./fec_outfiles/fec_out_");
	strcat(filename, record->sequence);
	if(record->conc != 0){
	  char* prot;
	  char* token;
      char* rest = record->input_filename;
	  int protargcount = 0;
      while ((token = strtok_r(rest, "/", &rest))){
		  prot = strtok_r(token, ".", &token);
	  }
	
      strcat(filename, conn);
	  strcat(filename, prot);
	  strcat(filename, conn);
	  strcat(filename, concname);
	}
	strcat(filename, end);

	FILE *out_file = fopen(filename, "w");
	
	for (force = (h/2.) + min_f; force <= max_f - (h/2.); force += delta_f)
	{
		index   = (int) ((force / delta_f)+0.5);
		h_ind   = (int) ((h / delta_f) + 0.5);
		exten   = -1*(energies[index + h_ind/2] - energies[index - h_ind/2]) / h;

		// print to output files 
		fprintf(out_file, "%f	", force); 
		fprintf(out_file, "%f	", exten);
		fprintf(out_file, "%f	\n", energies[index]);
	}
	
	fclose(out_file);
  
}
/* This version determines helix packing between all helices in a protein */

// dmf 6.27.17 
// ISSUE:  'axial' contact distances vulnerable to averaging over curved helices -
//         there may be work-arounds (e.g. introduce kinks in the .dssp file), or
//         there may be tests that could be added.
// ISSUES: hardwired location of ancillary files [need to work out a usage protocol]
// ISSUES: contact.txt list doesn't indicate which to which in the list. how to clarify? 
//         they are pairs, but tries to create unique residues list for the count.
//
// dmf 6.29.17
// FUTURE: would be easy to just to a brute-force screen of the distance
//         between the 'local' axes in order to find the pairs that are closest together.
//         once that point is found (which should generally agree with the more idealized
//         algorithm), create a short, localized axis (say over 2-turns of helix, or
//         whatever would be a typical packing area, so 11 residues or so (center +/- 5 residues)
//         and get the idealized crossing angles over that interface. SHOULD BE DONE. this would
//         actually give a more realistic picture of both packing distance and packing angle.
//
// SEE website for original code: http://fbs3pcu113.leeds.ac.uk/HelixPackingPair/source.html
// 
// creates the following files: 
// 	** helix_packing_pair.txt - simplest list of packing pairs 
// 	** helices.txt - lists helix numbers, packing angle and distance of contact pairs
// 	** contact.txt - lists residues of helix pairs involved in contacts 
// 	axis.txt - seems to be coordinates of local segment axes for each helix
// 	geom.txt - helix-fitting measures
// 	helix_shape.txt - lists shape of helices
//
// dmf 6.27.17 
// Fixes: 
// 1) change #define DLINLEN to 160 (was 150)
//  		Addresses: Problem arising from DSSP webserver line length issue, but OK now
// 2) change variable output_filename[20] to output_helices[24], or output_axis[24], or 
//    output_contact[24], output_geom[24] (in various functions), and substitute in code throughout;
//    also, output_file1[25] to output_packing[28]; also output_file2[20] to output_shape[24]
// 3) change usage of obpdbfile, so that it allows 1xyz.pdb as a default alternate to pdb1xyz.ent
// 4) added global flag new_open to control contact file opens (was flagged on helix # -> bug if
//    helix '0' reappeared in processing list)
// dmf 6.28.17
// 5) added output to axis.py with axial points and contact points for display in PyMol. see comments.
//          Addresses: Create a graphical translation of the tabular data by using axis.txt data
//          to do PyMol workaround like this: 
//				pseudoatom pt1, pos=[x1, y1, z1]
//				pseudoatom pt2, pos=[x2, y2, z2]
//				distance /pt1, /pt2
//				set dash_gap, 0
// 			see: https://sourceforge.net/p/pymol/mailman/message/25795427/
// dmf 6.29.17
// 6) changed the 'long helix' handling limit from 30 to 34
// 7) added segment_segment_closest_points3d() to get contact points within the helix segments;
//    added code to skew.h to do this; rearranged some code. see comments.
//    	   Addresses: Original close contact code determines helix distance based on infinite extension
//         of the helix axial lines. THIS IS REQUIRED FOR THE SIMPLE HELIX ANGLE
//         CALCULATION.  But this also leads to nonsensical 'contact points' and
//         likely is the source of the helix distances being too small. some helix distances 
//		   are maybe too small (6.7A, for example). Gly/Gly packing gives ~6.4A axial distance 
//		   (ballpark). 'close' packing may also be arising due to use of average helix axis 
// 		   over the whole (possibly curved) helix.
// dmf 7.12.17
// 8) recalculate helix_pair.distance so that the helix_packing_pair.txt output is consistent
//    with the numbers output to axis.py
// dmf 7.25.17 
// 9) begin working on adding a Identifier to the output file names.
//      1st - restructure main() code in anticipation of filename modification by pdb id
//      2nd - there are seven (7) output files to handle:
//              helix_shape.txt         - output_shape;     open in main() only
//              helix_packing_pair.txt  - output_packing;   open in main() only
//              helices.txt             - output_helices;   open in main() only
//              axis.py                 - pymol_axis;       open in main() and get_bending_angle()
//              axis.txt                - output_axis;      open in get_local_axis(), get_bending_angle() and two_helix_contact_vectors()
//              contact.txt             - output_contact;   open in atom_distance() only
//              geom.txt                - output_geom;      open in fit() only
//      3rd - added function to generate filenames with pdb_id appended - all good.
//          Addresses: output files are not named per the files being processed (e.g.
//          should have "1ffh_output.txt" or whatever, not just "output.txt")
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "skew.h"

// dmf 6.27.17
// #define DEBUG

// dmf 6.27.17 increased DLINLEN from 150
#define DLINLEN 160
#define MAXHELICES 1000               /* max helices in protein - set between 500-1000 */
#define MAXRESIDUES 200               /* maximum residues in a helix - set between 100-200 */
#define LINLEN 128
#define CLINLEN 450
#define MAXHELIXATOMS (MAXRESIDUES*10)
#define SQR(X) ((X)*(X))
#define TOLERANCE 0.6                 /* to be used in packing threshold calculation */
#define TOLERANCE2 (106.0/100.0)      /* to be used in determining atom-atom bonds   */
#define MIN_CONTACT_RESIDUES 3        /* minimum contacting residues per helix that constitute packing */
#define OBPDBDIR "/usr3/database/pdbobso/"

/* Global Variables */

// dmf 6.27.17
int new_open = 1;

// dmf 7.29.17
char output_helices[30];
char output_packing[30];
char output_shape[30];
char output_axis[30];
char output_geom[30];
char output_contact[30];
char pymol_axis[30];

/* For a helix containing 100 residues: there are 97 local axes, 98 local origins, 32 bending angles */

struct HELIX
{
   int helix_no;
   char pdb[5];
   char chain;    
   char residues[MAXRESIDUES+1];
   float residue_numbers[MAXRESIDUES];
   float ca_coord[MAXRESIDUES][3];
   double unit_local_axis[MAXRESIDUES-3][3];
   double origin[MAXRESIDUES-2][3];
   double bending_angle[MAXRESIDUES/3];
   double max_bending_angle;
   int residues_total;
   int atoms_total;
   char geometry;      /* S for short, U for unknown, K for kinked, C for curved, L for linear */
};

struct ATOM
{
   int atom_number;
   char atom_name[5];
   char residue_name[4];
   char chain;
   float residue_number;
   float atom_coord[3];
   int hdonor;                /* 1 for hbond donor or 0 if not */
   int charge;                /* 0 for neutral, 1 for positive, 2 for negative */
};

struct DISTANCE
{
   float distance;
   char atom_name1[5];
   char atom_name2[5];
   float residue_number1;
   float residue_number2;
   int helix_number1;
   int helix_number2;
};

struct HELIXPAIR
{
   int helix_one;
   int helix_two;
   int neighbours;            /* 0 for unrelated helices, 1 for neighbours */
   int packed;                /* 0 for unpacked, 1 for packed helices */
   int h1_residues;           /* number of residues of helix 1 involved in contact */
   int h2_residues;           /* number of residues of helix 2 involved in contact */
   int h1_start;              /* first residue of helix one in contact area, zero means first residue in helix (and not first in protein sequence) */
   int h2_start;              /* first residue of helix two in contact area, zero means first residue in helix */
   int h1_end;                /* last residue of helix one contact area */
   int h2_end;                /* last residue of helix two contact area */
   int vdw;                   /* number of interhelical vdw contacts */
   int hbond;                 /* number of interhelical hbond contacts */
   int electrostatic;         /* number of interhelical electrostatic interactions */
   int covalent;              /* number of interhelical covalent bonds */     
   double angle1;             /* the interhelical dihedral angle using all helical axis vectors */
   double angle2;             /* the interhelical dihedral angle using only the axis vectors in the contact area */
   double distance;           /* length of line of closest approach between the 2 helix axes using first method */
};
   
/* Prototypes */

struct HELIX* read_helices(FILE *fpi_dssp, int *helices_total, char *pdb_id);
struct ATOM** read_atom(FILE *fpi_pdb, struct HELIX*, int *helices_total, int *helices_atom_total);
void get_atom_info(struct ATOM**, struct HELIX*, int *helices_total);
void get_ca_coords(struct HELIX*, struct ATOM**, int *helices_total);
void get_local_axis(int helix_number, struct HELIX*);
void get_bending_angle(int helix_number, struct HELIX*);
void fit(int helix_number, struct HELIX*);
double** matinv3(double **h);
double** matinv2(double **p);
void destroy_matrix3(double **s);
void destroy_matrix2(double **q);
struct HELIXPAIR** neighbours(int *helices_total); 
struct DISTANCE** residue_distance(int helix1, int helix2, struct HELIX*, struct HELIXPAIR**);
struct DISTANCE** atom_distance(int helix1, int helix2, struct HELIX*, struct ATOM**, struct HELIXPAIR**);
void destroy_residue_residue(struct DISTANCE**, struct HELIX*, int helix_number);
void destroy_atom_atom(struct DISTANCE**, struct HELIX*, int helix_number);
void two_helix_all_vectors(int helix_A, int helix_B, struct HELIX*, struct HELIXPAIR**);
void two_helix_contact_vectors(int helix_A, int helix_B, struct HELIX*, struct HELIXPAIR**);
void destroy_helix_atom(struct ATOM**, int *helices_total);
void destroy_helix_pair(struct HELIXPAIR**, int *helices_total);
// dmf 7.29.17
void create_filenames(char *pdb_id);

/* --------------------------------Entry Point---------------------------- */

int main(void)     
{
   /* Variables */

   FILE *fpi_dssp;
   FILE *fpi_pdb;
   FILE *fpi_cath;
   FILE *fpo_helices;
   FILE *fpo_packing;
   FILE *fpo_shape;
   struct HELIX *helix;
   struct ATOM **helix_atom;
   struct DISTANCE **residue_residue;
   struct HELIXPAIR **helix_pair;
   struct DISTANCE **atom_atom;
   int g,h,i,j;
   int helices_total;
   int helices_atom_total;
   char pdb_id[5];
   char temp_pdb[5]="";
   char dsspfile[50];
   char PDBDIR[39];
   char DSSPDIR[41];
   char pdbfile[50];
   char obpdbfile[40];
   char cathfile[50]="";
   char line[CLINLEN];
    
    // dmf 6.30.17
    FILE *fpo_pyaxis;

#ifdef DEBUG
	printf("\n\tHello world!\n\tDebugging mode active\n\n"); 
#endif

   printf("\nInput filename read by taking first four characters of each line.\n");
   printf("Four characters: <pdb code>\n\n");

   printf("Please enter PDB directory (or leave blank for current directory): ");
   fgets (PDBDIR,39,stdin);

   if(PDBDIR[strlen(PDBDIR)-1]=='\n') PDBDIR[strlen(PDBDIR)-1]='\0';
   if(strlen(PDBDIR) && PDBDIR[strlen(PDBDIR)-1]!='/') strcat(PDBDIR,"/");
   
   printf("Please enter DSSP directory (or leave blank for current directory): ");
   fgets (DSSPDIR,41,stdin);

   if(DSSPDIR[strlen(DSSPDIR)-1]=='\n') DSSPDIR[strlen(DSSPDIR)-1]='\0';
   if(strlen(DSSPDIR) && DSSPDIR[strlen(DSSPDIR)-1]!='/') strcat(DSSPDIR,"/");

   printf("Please enter the input filename (eg cath.txt, Hreps.v2.4, Nreps.v2.4): ");
//   scanf ("%s",cathfile);
// dmf 
	sprintf(cathfile,"test.txt");   
 
//   while(getchar() != 10); /* Clear the input stream */
   printf("\n");

   if((fpi_cath=fopen(cathfile,"r"))==NULL)
   {
      printf("Error opening %s\n",cathfile);
      exit(1);
   }

// ***** cath.txt read loop *****
// the following reads in each line from the cathfile input list and processes it
   while(!feof(fpi_cath)) 
   {
      fgets(line, CLINLEN, fpi_cath);
       
      if(line[0]=='\n' || line[0]==' ') continue;
      if(line[strlen(line)-1]=='\n') line[strlen(line)-1]='\0';

// dmf 7.25.17 want to grab the 'line' to a filename specifier to be 
// added to each of the output filenames.
// dmf 7.25.17 this grabs the line from the input list
      strncpy(pdb_id,line,4);
// dmf 7.25.17 this truncates the line to 4 characters + \0
      pdb_id[4]='\0';

// dmf 7.25.17 this compares the new name to the current name. Note how it works:
//  The strcmp() and strncmp() functions return an integer greater than,
//  equal to, or less than 0, according as the string s1 is greater than,
//	equal to, or less than the string s2. 
      if(strcmp(temp_pdb,pdb_id))
      {
         strcpy(temp_pdb,pdb_id);

         strcpy(dsspfile,DSSPDIR);
         strcat(dsspfile,pdb_id);
         strcat(dsspfile,".dssp");

         strcpy(pdbfile,PDBDIR);
         strcat(pdbfile,"pdb");
         strcat(pdbfile,pdb_id);
         strcat(pdbfile,".ent");

// dmf 6.28.17 change to allow 1xyz.pdb as a default alternate to pdb1xyz.ent
         strcpy(obpdbfile,PDBDIR);
         strcat(obpdbfile,pdb_id);
         strcat(obpdbfile,".pdb");

         printf("\nAnalysis of %s in progress\n",pdb_id);
         printf("Input Files: %s, %s\n\n",dsspfile,pdbfile);

          // dmf 7.29.17 create the output filenames
          create_filenames(pdb_id);
          
// ** moved from above **
          
          // dmf 7.25.17 - want to modify output_packing to include identifying string.
          // therefore, this statement needs to be moved to after file input, below.
          if((fpo_packing=fopen(output_packing, "w"))==NULL)
          {
              printf("\n\n** Error writing to file '%s'!",output_packing);
              exit(1);
          }
          
          fprintf(fpo_packing,"Protein\tHelix1\tHelix2\tCont 1\tCont 2\tGlobal Angle\tLocal Angle\tDistance\tCovalnt\tElectro\tH-Bond\tVDW\n");
          
          // dmf 7.25.17 - want to modify output_shape to include identifying string.
          // therefore, this statement needs to be moved to after file input, below.
          if((fpo_shape=fopen(output_shape, "w"))==NULL)
          {
              printf("\n\n** Error writing to file '%s'!",output_shape);
              exit(1);
          }
          
          fprintf(fpo_shape,"Protein\tChain\tHelix\tLength\tGeom\tMax Bending Angle\n");
          
// ** end of moved from above **

// dmf 7.25.17 - want to modify output_helices to include identifying string. 
         if((fpo_helices=fopen(output_helices, "w"))==NULL)
         {
            printf("\n\n** Error writing to file '%s'!",output_helices);
            exit(1);
         }

         if((fpi_dssp=fopen(dsspfile,"r"))==NULL)
         {
            printf("\n\nError opening %s\n",dsspfile);
            exit(1);
         }

         helix=read_helices(fpi_dssp, &helices_total, pdb_id);

         fclose(fpi_dssp);

         if((fpi_pdb=fopen(pdbfile,"r"))==NULL)
         {
            if((fpi_pdb=fopen(obpdbfile,"r"))==NULL)
            {
               printf("\n\nError opening %s\n",pdbfile);
               printf("Error opening %s\n",obpdbfile);
               exit(1);
            }
         }

         helix_atom=read_atom(fpi_pdb, helix, &helices_total, &helices_atom_total);

         fclose(fpi_pdb);

         get_atom_info(helix_atom, helix, &helices_total);

         get_ca_coords(helix, helix_atom, &helices_total);

         for(i=0;i<helices_total;i++)
         {
            get_local_axis(i, helix);

            get_bending_angle(i, helix);
       
            fit(i, helix);
         } 
   
         for(i=0;i<helices_total;i++)
         {
            fprintf(fpo_helices,"protein: %s, chain: %c, helix: %d, start residue: %.2f, last residue: %.2f\n",helix[i].pdb,helix[i].chain,helix[i].helix_no,helix[i].residue_numbers[0],helix[i].residue_numbers[helix[i].residues_total-1]);
            fprintf(fpo_helices,"number of residues: %d, sequence: %s\n",helix[i].residues_total,helix[i].residues);
            fprintf(fpo_helices,"maximum bending angle: %f degrees, overall geometry: %c\n\n",helix[i].max_bending_angle,helix[i].geometry);
         }

         fprintf(fpo_helices,"Total Number of Helices = %d\n\n",helices_total);

         helix_pair=neighbours(&helices_total);

         fprintf(fpo_helices,"Neighbouring Helices\n\n");

         for(j=0;j<helices_total;j++)
         {
            printf(".");
            fflush(stdout);
    
            for(i=0;i<=j;i++)
            {
               residue_residue=residue_distance(i, j, helix, helix_pair);
          
               destroy_residue_residue(residue_residue, helix, i);

               if(helix_pair[i][j].neighbours==1)
               {
                  atom_atom=atom_distance(i, j, helix, helix_atom, helix_pair);

                  fprintf(fpo_helices,"helix %d & ",helix_pair[i][j].helix_one);
                  fprintf(fpo_helices,"helix %d  ",helix_pair[i][j].helix_two);
                  fprintf(fpo_helices,"neighbours: %d  ",helix_pair[i][j].neighbours);
                  fprintf(fpo_helices,"packed: %d  ",helix_pair[i][j].packed);
                  fprintf(fpo_helices,"helix %d contact residues: %d  ",helix_pair[i][j].helix_one,helix_pair[i][j].h1_residues);
                  fprintf(fpo_helices,"helix %d contact residues: %d  ",helix_pair[i][j].helix_two,helix_pair[i][j].h2_residues);
                  fprintf(fpo_helices,"helix %d first contact residue: %d  ",helix_pair[i][j].helix_one,helix_pair[i][j].h1_start);
                  fprintf(fpo_helices,"last contact residue: %d  ",helix_pair[i][j].h1_end);
                  fprintf(fpo_helices,"helix %d first contact residue: %d  ",helix_pair[i][j].helix_two,helix_pair[i][j].h2_start);
                  fprintf(fpo_helices,"last contact residue: %d  ",helix_pair[i][j].h2_end);
                  fprintf(fpo_helices,"covalents: %d  ",helix_pair[i][j].covalent);
                  fprintf(fpo_helices,"electrostatics: %d  ",helix_pair[i][j].electrostatic);
                  fprintf(fpo_helices,"hbonds: %d  ",helix_pair[i][j].hbond);
                  fprintf(fpo_helices,"vdws: %d\n",helix_pair[i][j].vdw);
                  destroy_atom_atom(atom_atom, helix, i);
               }
            }
         }

         fprintf(fpo_helices,"\nPacked Helices: Angles & Distance of Closest Approach\n\n");

         for(j=0;j<helices_total;j++)
         {
            for(i=0;i<=j;i++)             
            {
               if((helix_pair[i][j].packed==1) && (helix[i].residues_total>=4) && (helix[j].residues_total>=4))
               {
                  two_helix_all_vectors(i, j, helix, helix_pair);
             
                  two_helix_contact_vectors(i, j, helix, helix_pair);

                  fprintf(fpo_helices,"Helix %d & Helix %d\n",i,j);
                  fprintf(fpo_helices,"Global Angle (from all vectors): %f degrees\nLocal Angle (from contact vectors): %f degrees\nInteraxial Distance: %f Angstroms\n\n",helix_pair[i][j].angle1,helix_pair[i][j].angle2,helix_pair[i][j].distance);
               }
            }
         }

         fprintf(fpo_helices,"\nAtomic List\n\n");

         for(i=0;i<helices_total;i++)
         {
            for(j=0;j<helix[i].atoms_total;j++)
            {
               fprintf(fpo_helices,"atom: %d,%s  hbond donor: %d  charge: %d  residue: %s  residue number: %.2f  chain: %c\n",helix_atom[i][j].atom_number,helix_atom[i][j].atom_name,helix_atom[i][j].hdonor,helix_atom[i][j].charge,helix_atom[i][j].residue_name,helix_atom[i][j].residue_number,helix_atom[i][j].chain);
            }

            fprintf(fpo_helices,"\n");
         }

         fprintf(fpo_helices,"total number of helical atoms = %d\n\n",helices_atom_total);

         fclose(fpo_helices);

         for(j=0;j<helices_total;j++)
         {
            for(i=0;i<=j;i++)             
            {
               if((helix_pair[i][j].packed==1) && (helix[i].residues_total>=4) && (helix[j].residues_total>=4))
               {
                  // output to "output_packing", file is already open with header text
                  fprintf(fpo_packing,"%s\t%d\t%d\t%d\t%d\t",helix[i].pdb,helix[i].helix_no,helix[j].helix_no,helix_pair[i][j].h1_residues,helix_pair[i][j].h2_residues); 
                  fprintf(fpo_packing,"%f\t%f\t%f\t",helix_pair[i][j].angle1,helix_pair[i][j].angle2,helix_pair[i][j].distance);
                  fprintf(fpo_packing,"%d\t%d\t%d\t%d\n",helix_pair[i][j].covalent,helix_pair[i][j].electrostatic,helix_pair[i][j].hbond,helix_pair[i][j].vdw);
               }
            }      
         }

         fclose(fpo_packing);

         for(i=0;i<helices_total;i++)
         {
            // output to "output_shape", file is already open with header text
            fprintf(fpo_shape,"%s\t%c\t%d\t%d\t",helix[i].pdb,helix[i].chain,helix[i].helix_no,helix[i].residues_total);
            fprintf(fpo_shape,"%c\t%f\n",helix[i].geometry,helix[i].max_bending_angle);     
         }      

         fclose(fpo_shape);

         printf("  Done\n");

         free(helix);

         destroy_helix_atom(helix_atom, &helices_total);

         destroy_helix_pair(helix_pair, &helices_total);
          
         // add tail information to axis.py
         // dmf 7.25.17 - want to modify pymol_axis to include identifying string
         // dmf 7.27.17 note that this file is opened and appended to in several different
         // functions, so this "a+" is required here.
         if((fpo_pyaxis=fopen(pymol_axis, "a+"))==NULL)
         {
            printf("\n\n** Error appending to file '%s'!",pymol_axis);
            exit(1);
         }
         fprintf(fpo_pyaxis,"set dash_gap, 0, cont*\n");
         fprintf(fpo_pyaxis,"set dash_radius, 0.40\n");
         fprintf(fpo_pyaxis,"set dash_round_ends, 0\n");
         fprintf(fpo_pyaxis,"set dash_color, 0xffcc00, dist*\n");
         fprintf(fpo_pyaxis,"hide labels, dist*\n");
         fclose(fpo_pyaxis);

      }
   }
// end of cathfile input read and process loop
// ***** end of cath.txt read loop *****
    
   fclose(fpi_cath);

   return 0;
}
// end of main();

/* ------------------------------------------------------------------------- */
      
/* Function to read DSSP file, get residues in helices, and initialise helix array */
struct HELIX* read_helices(FILE *fpi_dssp, int *helices_total, char *pdb_id)
{
   /* Variables */

   char line[DLINLEN];
   struct HELIX *helix;
   char res[7];
   char res_sub_type;
   float res_number;
   float res_sub_number;
   int g,h,i=0,j,k=0;
   char previous_structure='Z';
   char current_structure='X';


   helix=(struct HELIX *) calloc(MAXHELICES,sizeof(struct HELIX));

   /* fill the array of the structure HELIX with junk for debugging */ 

   for (g=0; g<MAXHELICES; g++)
   {
      helix[g].residues_total=-1;
      helix[g].helix_no=-1;
      helix[g].atoms_total=-1;
      helix[g].max_bending_angle=0;
      helix[g].chain='Z';
      helix[g].geometry='S';

      for (h=0; h<MAXRESIDUES; h++)
      {
         helix[g].residues[h]='Z';
         helix[g].residue_numbers[h]=-1;
         helix[g].ca_coord[h][0]=-9999;
         helix[g].ca_coord[h][1]=-9999;
         helix[g].ca_coord[h][2]=-9999;
      }
      helix[g].residues[h+1]='\0';

      for (h=0; h<MAXRESIDUES-3; h++)
      {
         helix[g].unit_local_axis[h][0]=-1;    
         helix[g].unit_local_axis[h][1]=-1;    
         helix[g].unit_local_axis[h][2]=-1;
      }

      for (h=0; h<MAXRESIDUES-2; h++)
      {
         helix[g].origin[h][0]=-1;
         helix[g].origin[h][1]=-1;
         helix[g].origin[h][2]=-1; 
      }

      for (h=0; h<MAXRESIDUES/3; h++)
      {
         helix[g].bending_angle[h]=-1;
      }         
   }

   /* start getting the dssp file line by line, and get to important bit */

   do
   {
      fgets(line, DLINLEN, fpi_dssp); 
   }
   while(!feof(fpi_dssp) && line[2]!='#');
  
   /* i is number of helices */
   /* k is number of residues in each helix */

   while(!feof(fpi_dssp))
   {
      fgets(line, DLINLEN, fpi_dssp);
  
      if(line[strlen(line)-1]=='\n') line[strlen(line)-1]='\0';

      if(line[13]!='!')
      {
         if(line[11]==' ') line[11]='0';
 
         for(j=0; j<6; j++) res[j]=line[j+5];  
         res[6]='\0';
        
         /* get residue number, e.g. residue 10 becomes 10.00 */

         res_number=atof(res);
       
         res_sub_type=' ';
         res_sub_number=0.0;

         /* get residue sub-label if it exists */

         if(line[10]!=' ' && line[10]!='0' && line[10]!='1' && line[10]!='2' && line[10]!='3' && line[10]!='4' && line[10]!='5' && line[10]!='6' && line[10]!='7' && line[10]!='8' && line[10]!='9')
         {
            res_sub_type=line[10];
            res_sub_number=res_sub_type-64;
            res_sub_number=res_sub_number/100;
         }

         /* this has converted residue 1A into 1.01, residue 2B into 2.02 etc. */

         res_number+=res_sub_number;

         current_structure=line[16];    /* update current secondary structure type */


         if(current_structure=='I' || current_structure=='G' || current_structure=='H')    /* if secondary structure of this line is a helix... */
         {
            strcpy(helix[i].pdb,pdb_id);

            helix[i].chain=line[11];
                  
            helix[i].residue_numbers[k]=res_number;
         
            helix[i].residues[k]=line[13];
            helix[i].residues[k+1]='\0'; 

            /* and increment k (the helix residue number) for the next residue of existing helix */

            k++;
       
            if(k==MAXRESIDUES)
            {
               printf("\n\nMaximum helix residue number has been reached for helix %d\n",i);
               exit(1);
            }
         }

         /* if secondary structure of this line is not a helix and the structure from */
         /* the previous line was... then assume end of helix and reset k (helix residue number) */
         /* and increment i (the helix number) in preparation for start of next helix. */

         else if((current_structure!='I' || current_structure!='G' || current_structure!='H') && (previous_structure=='I' || previous_structure=='G' || previous_structure=='H'))
         {
            helix[i].residues_total=k;

            helix[i].helix_no=i;

            k=0;

            i++;        

            if(i==MAXHELICES)
            {
               printf("\n\nMaximum number of helices allowed has been reached\n");
               exit(1);
            }
         }

         previous_structure=current_structure;         
      }

      else if((line[13]=='!') && (previous_structure=='I' || previous_structure=='G' || previous_structure=='H'))
      {
         previous_structure='Z';
         current_structure='X';

         helix[i].residues_total=k;

         helix[i].helix_no=i;

         k=0;

         i++;

         if(i==MAXHELICES)
         {
            printf("\n\nMaximum number of helices allowed has been reached\n");
            exit(1);
         }
      }                     
   }

   *helices_total=i;

   return helix;
}

/* ------------------------------------------------------------------------- */

/* Function to read PDB file and get atom details if they are in the DSSP defined helices */
struct ATOM** read_atom(FILE *fpi_pdb, struct HELIX *helix, int *helices_total, int *helices_atom_total)
{
   /* Variables */

   struct ATOM **helix_atom;
   char line[LINLEN];
//   char coord[9];
// dmf 6.27.17
   char coord[10];
   char resseq[6];
   char resname[4];
   char chain;
   char res_sub_type;
   float res_sub_number;
   float current_residue_number=-1.5;
   float last_residue_in_helix=-1.5;
   char number[6];
   int g,h,i=0,j=0,k,l,m=0,n=0;


   /* Allocate the memory for the 2d array dynamically
   *
   * -- max atoms per helix (1000-2000)-->
   *        |
   *  number of helices
   *        \/
   */
 
   helix_atom = (struct ATOM **) calloc(*helices_total,sizeof(struct ATOM*));
   
   for(k=0; k<*helices_total; k++)
   {
      helix_atom[k] = (struct ATOM *) calloc(MAXHELIXATOMS,sizeof(struct ATOM));
   }

   /* fill all atom and residue numbers in entire 2-D array with junk */

   for(g=0; g<*helices_total; g++)
   {
      for(h=0; h<MAXHELIXATOMS; h++)
      {            
         helix_atom[g][h].atom_number=-1;
         helix_atom[g][h].residue_number=-1;
         helix_atom[g][h].chain='Z';
         helix_atom[g][h].atom_coord[0]=-9999;
         helix_atom[g][h].atom_coord[1]=-9999;
         helix_atom[g][h].atom_coord[2]=-9999;
         helix_atom[g][h].hdonor=0;
         helix_atom[g][h].charge=0;
      }   
   }

   /* Read input file until all helices are completed */

   while(!feof(fpi_pdb) && strncmp(line,"END",3))
   {
      fgets(line, LINLEN, fpi_pdb);

      if(line[strlen(line)-1]=='\n') line[strlen(line)-1]='\0';

      /* If record name is ATOM or HETATM and the atom name is not H ... */

      if((line[13]!='H') && (!strncmp(line,"ATOM  ",6) || !strncmp(line,"HETATM",6)))
      {
         if(line[21]==' ') line[21]='0';
         chain=line[21];

         for(k=0;k<3;k++) resname[k]=line[k+17];	
         resname[3]='\0';

         /* if the chain is right and the residue is not a water molecule */

         if((chain==helix[i].chain) && (strcmp(resname,"HOH")))        
         {
            /* current residue number from PDB atom list */
            /* e.g. residue 1 becomes 1.00 */

            for(k=0;k<5;k++) resseq[k]=line[k+22];
	
            resseq[5]='\0';
            current_residue_number=atof(resseq);

            res_sub_type=' ';
            res_sub_number=0.0;

            /* get residue sub-label if it exists */

            if(line[26]!=' ' && line[26]!='0' && line[26]!='1' && line[26]!='2' && line[26]!='3' && line[26]!='4' && line[26]!='5' && line[26]!='6' && line[26]!='7' && line[26]!='8' && line[26]!='9') 
            {
               res_sub_type=line[26];
               res_sub_number=res_sub_type-64;
               res_sub_number=res_sub_number/100;
            }

            /* this has converted residue 1A into 1.01, residue 1Z into 1.26 etc. */

            current_residue_number+=res_sub_number;

            /* For within a helix, where j (helix residue number) must be incremented before assignments to atom variables are made */

            if((current_residue_number==helix[i].residue_numbers[j+1]) && (last_residue_in_helix==helix[i].residue_numbers[j]))
            {
               j++;
            }

            /* For start or within helix, where PDB residue number equals helix residue number */

            if(current_residue_number==helix[i].residue_numbers[j])
            {
               /* atom_number */

               for(k=0;k<5;k++) number[k]=line[k+6];
    
               number[5]='\0';
               helix_atom[i][m].atom_number=atoi(number);
	
               /* atom name */

               for(k=0;k<4;k++) helix_atom[i][m].atom_name[k]=line[k+12];
	
               helix_atom[i][m].atom_name[4]='\0';

               /* residue name */
        
               strcpy(helix_atom[i][m].residue_name,resname);

               /* chain */

               helix_atom[i][m].chain=chain;
                  	
               /* residue number */

               helix_atom[i][m].residue_number=current_residue_number;

               /* co-ordinates */

               k=30;
              
               for(g=0;g<3;g++)
               {

                  for(h=0;h<9;h++) coord[h]=line[h+k];
               
                  coord[9]='\0';
               
                  helix_atom[i][m].atom_coord[g]=atof(coord);

                  k=k+8;
               }

               last_residue_in_helix=current_residue_number; /* update last residue number in helix */

               m++;    /* m represents total atoms in a helix */

               n++;    /* n represents total atoms in all helices */
            }

            /* for a search for the start of another helix */ 

            if((current_residue_number!=helix[i].residue_numbers[j]) && (last_residue_in_helix==helix[i].residue_numbers[j]))
            {
               j++;
            }

            /* when a helix is at an end */

            if(j==helix[i].residues_total)
            {
               helix[i].atoms_total=m;
               i++;
               j=0;  
               m=0; 
            }

            /* when all helices are completed */

            if(i==*helices_total)
            {
               break;
            }            
         }
      }
   }

  *helices_atom_total=n;        

  return helix_atom;
}

/* ------------------------------------------------------------------------- */

/* Function to update ATOM structures in the helix_atom 2-D array with hydrogen bond donor and electrostatic info */
void get_atom_info(struct ATOM **helix_atom, struct HELIX *helix, int *helices_total)
{
   /* Variables */

   int i,j,k,l;
   FILE *fp;            // dmf 7.27.17 - leave as fp
   char line[LINLEN];
   char residue[4];
   char atom[5];
   char hbond[2];
   char electro[2];
   int hbond_donor;
   int electrostatic;


   if((fp=fopen("translation.txt","r")) == NULL)
   {
      printf("\n\nError opening translation.txt\n");
      exit(1);
   }
  
   while(!feof(fp))
   {
      fgets(line,LINLEN,fp);

      if(line[strlen(line)-1]=='\n') line[strlen(line)-1]='\0';
    
      for(k=0;k<3;k++) residue[k]=line[k];
      residue[3]='\0';
    
      for(l=0;l<4;l++) atom[l]=line[4+l];     
      atom[4]='\0';
    
      hbond[0]=line[9];                      
      hbond[1]='\0';
      hbond_donor=atoi(hbond);
       
      electro[0]=line[11];                   
      electro[1]='\0';
      electrostatic=atoi(electro);
     
      for(i=0;i<*helices_total;i++)
      {
         for(j=0;j<helix[i].atoms_total;j++)
         {
            if((!strcmp(residue,helix_atom[i][j].residue_name)) && (!strcmp(atom,helix_atom[i][j].atom_name)))
            {
               helix_atom[i][j].hdonor=hbond_donor;
               helix_atom[i][j].charge=electrostatic;
            }
         }
      }
   }
// printf("\n\nresidue %s, atom %s, hbond %d, charge %d\n",residue,atom,hbond_donor,electrostatic);

   fclose(fp);
}

/* ------------------------------------------------------------------------- */

/* Little function to get CA co-ordinates from one structure to another */
void get_ca_coords(struct HELIX *helix, struct ATOM **helix_atom, int *helices_total)
{  
   /* Variables */

   int i,j,k;
   float current_residue;
   float previous_residue;


   for (i=0;i<*helices_total;i++)
   {
      k=0;              /* set back to zero for new helix */

      previous_residue=-1.5;   /* reset for new helix */

      for(j=0;j<helix[i].atoms_total;j++)
      {
         current_residue=helix_atom[i][j].residue_number;

         if((!strcmp(helix_atom[i][j].atom_name," CA ")) && (current_residue!=previous_residue))
         {
            helix[i].ca_coord[k][0]=helix_atom[i][j].atom_coord[0];
            helix[i].ca_coord[k][1]=helix_atom[i][j].atom_coord[1];
            helix[i].ca_coord[k][2]=helix_atom[i][j].atom_coord[2];

            k++;        /* k is number of CA atoms */

            previous_residue=current_residue;
         }

         if(k==helix[i].residues_total+1)
         {
            printf("\n\n*** Maximum number of helix CA co-ordinates transferred ***\n");
            exit(1);
         }
      }
   }
}

/* ------------------------------------------------------------------------- */

/* Function to get the unit local axes (upto 97 in a 100 residue helix) of a single helix */
/* and then to get the local helix origins (upto 98 in a 100 residue helix */
void get_local_axis(int helix_number, struct HELIX *helix)
{  
   /* Variables */

   int i,j;
   double vec12[3];
   double vec23[3];
   double vec34[3];
   double dv13[3];
   double dv24[3];
   double cross_product[3];
   double dot_product;
   double mag_cross_product;
   double dmag;
   double emag;
   double costheta;
   double costheta1;
   double radmag;
   double rad[3];
   FILE *fpo_axis;

   i=helix_number;

   if(i==0)
   {
// dmf 7.25.17 - want to modify output_axis to include identifying string. 
      if((fpo_axis=fopen(output_axis, "w"))==NULL)
      {
         printf("\n\n** Error writing to file '%s'!",output_axis);
         exit(1);
      }
   }
   else   
   {
// dmf 7.25.17 - want to modify output_axis to include identifying string. 
      if((fpo_axis=fopen(output_axis, "a+"))==NULL)
      {
         printf("\n\n** Error appending to file '%s'!",output_axis);
         exit(1);
      }
   }
    
   fprintf(fpo_axis,"Helix Number %d\n\n",i);

// dmf 6.27.17 need to start assembling the PyMol script file here...
    
   if(helix[i].residues_total>=4)
   {
      for(j=0;j<helix[i].residues_total-3;j++)
      {
         /* get 4 consecutive CA atoms */
         /* j = CA atom number */
         /* [0] = X-axis, [1] = Y-axis, [2] = Z-axis */
         /* vectors joining CA atoms */

         if(helix[i].ca_coord[j][0]==-9999 || helix[i].ca_coord[j][1]==-9999 || helix[i].ca_coord[j][2]==-9999 || helix[i].ca_coord[j+1][0]==-9999 || helix[i].ca_coord[j+1][1]==-9999 || helix[i].ca_coord[j+1][2]==-9999 || helix[i].ca_coord[j+2][0]==-9999 || helix[i].ca_coord[j+2][1]==-9999 || helix[i].ca_coord[j+2][2]==-9999 || helix[i].ca_coord[j+3][0]==-9999 || helix[i].ca_coord[j+3][1]==-9999 || helix[i].ca_coord[j+3][2]==-9999)
         {
            printf("\n\n** error - c-alpha atom limit breached in vector analysis\n");
            exit(1);
         }

         vec12[0]=helix[i].ca_coord[j+1][0] - helix[i].ca_coord[j][0];
         vec12[1]=helix[i].ca_coord[j+1][1] - helix[i].ca_coord[j][1];
         vec12[2]=helix[i].ca_coord[j+1][2] - helix[i].ca_coord[j][2];

         vec23[0]=helix[i].ca_coord[j+2][0] - helix[i].ca_coord[j+1][0];
         vec23[1]=helix[i].ca_coord[j+2][1] - helix[i].ca_coord[j+1][1];
         vec23[2]=helix[i].ca_coord[j+2][2] - helix[i].ca_coord[j+1][2];

         vec34[0]=helix[i].ca_coord[j+3][0] - helix[i].ca_coord[j+2][0];
         vec34[1]=helix[i].ca_coord[j+3][1] - helix[i].ca_coord[j+2][1];
         vec34[2]=helix[i].ca_coord[j+3][2] - helix[i].ca_coord[j+2][2];

         /* difference of vectors joining CA atoms */
         /* these vectors are perpendicular to the helix axis */

         dv13[0]=vec12[0] - vec23[0];
         dv13[1]=vec12[1] - vec23[1];
         dv13[2]=vec12[2] - vec23[2];

         dv24[0]=vec23[0] - vec34[0];
         dv24[1]=vec23[1] - vec34[1];
         dv24[2]=vec23[2] - vec34[2];

         /* x, y, z components of the cross-product vector (the axis) */

         cross_product[0]=(dv13[1] * dv24[2]) - (dv13[2] * dv24[1]);
         cross_product[1]=(dv13[2] * dv24[0]) - (dv13[0] * dv24[2]);
         cross_product[2]=(dv13[0] * dv24[1]) - (dv13[1] * dv24[0]);

         /* get unit vectors of the axis */

         mag_cross_product=sqrt(SQR(cross_product[0])+SQR(cross_product[1])+SQR(cross_product[2]));

         cross_product[0]=cross_product[0]/mag_cross_product;
         cross_product[1]=cross_product[1]/mag_cross_product;
         cross_product[2]=cross_product[2]/mag_cross_product;
         
         helix[i].unit_local_axis[j][0]=cross_product[0];
         helix[i].unit_local_axis[j][1]=cross_product[1];
         helix[i].unit_local_axis[j][2]=cross_product[2];

         fprintf(fpo_axis,"unit local axis %d: X %f, Y %f, Z %f\n",j,helix[i].unit_local_axis[j][0],helix[i].unit_local_axis[j][1],helix[i].unit_local_axis[j][2]);
     
         dmag=sqrt(SQR(dv13[0]) + SQR(dv13[1]) + SQR(dv13[2]));
         emag=sqrt(SQR(dv24[0]) + SQR(dv24[1]) + SQR(dv24[2]));
         dot_product=(dv13[0]*dv24[0]) + (dv13[1]*dv24[1]) + (dv13[2]*dv24[2]);
         costheta=dot_product/(dmag*emag);

         /* radius of local helix cylinder */

         costheta1=1.0-costheta;
         radmag=sqrt(dmag*emag)/(costheta1*2.0);

         /* local helix origins */
         /* convert dv13 and dv24 to unit vectors */

         dv13[0]=dv13[0]/dmag;
         dv13[1]=dv13[1]/dmag;
         dv13[2]=dv13[2]/dmag;

         /* convert vectors to length of cylinder radius */

         rad[0]=radmag*dv13[0];
         rad[1]=radmag*dv13[1];
         rad[2]=radmag*dv13[2];

         /* get two helix origins for each local helix axis */
         /* this is the first.... */

         helix[i].origin[j][0]=helix[i].ca_coord[j+1][0]-rad[0];
         helix[i].origin[j][1]=helix[i].ca_coord[j+1][1]-rad[1];
         helix[i].origin[j][2]=helix[i].ca_coord[j+1][2]-rad[2];

         dv24[0]=dv24[0]/emag;
         dv24[1]=dv24[1]/emag;
         dv24[2]=dv24[2]/emag;

         rad[0]=radmag*dv24[0];
         rad[1]=radmag*dv24[1];
         rad[2]=radmag*dv24[2];

         /* and the second....  */
         /* this is actually recorded over by the first origin of the next local axis, except at the very last axis! */

         helix[i].origin[j+1][0]=helix[i].ca_coord[j+2][0]-rad[0];
         helix[i].origin[j+1][1]=helix[i].ca_coord[j+2][1]-rad[1];
         helix[i].origin[j+1][2]=helix[i].ca_coord[j+2][2]-rad[2];

         fprintf(fpo_axis,"helix origin %d: X %f, Y %f, Z %f\n",j,helix[i].origin[j][0],helix[i].origin[j][1],helix[i].origin[j][2]);
         fprintf(fpo_axis,"helix origin %d: X %f, Y %f, Z %f\n\n",j+1,helix[i].origin[j+1][0],helix[i].origin[j+1][1],helix[i].origin[j+1][2]);
          

      }
   }
   else
   {
      fprintf(fpo_axis,"Helix %d is less than 4 residues and has no axis\n\n",i);
   }

   fclose(fpo_axis);
}

/* ------------------------------------------------------------------------- */

/* Function to get the bending angle between two local axes of a single helix (angles between axes j--j+3, j+3--j+6, etc.) and then get the maximum bending angle in the helix */
void get_bending_angle(int helix_number, struct HELIX *helix)
{
   /* Variables */

   int i,j,k=0,l;
   FILE *fpo_axis;
   double angle;
   double max_angle=0;
   double pi;

   // dmf 6.27.17 will add a axis.py file that will contain a workaround
   // to display the axis elements for each helix in PyMol. Will also want
   // to add a vector for the point of close approach.
   FILE *fpo_pyaxis;
   int vec_flag;
   char pt_label[20], previous_pt_label[20];
   char dLabel[20]; 

   i=helix_number;

   pi=180.0/acos(-1.0);

// dmf 7.25.17 - want to modify output_axis to include identifying string.
   if((fpo_axis=fopen(output_axis, "a+"))==NULL)
   {
      // dmf 6.28.17 note that this is the same files as output_axis, so exists
      printf("\n\n** Error appending to file '%s'!",output_axis);
      exit(1);
   }

    if(i==0)
    // dmf 6.28.17 handle initialization and appending to pymol output file
    {
// dmf 7.25.17 - want to modify pymol_axis to include identifying string. 
        if((fpo_pyaxis=fopen(pymol_axis, "w"))==NULL)
        {
            printf("\n\n** Error writing to file '%s'!",pymol_axis);
            exit(1);
        }
    }
    else
    {
// dmf 7.25.17 - want to modify pymol_axis to include identifying string. 
        if((fpo_pyaxis=fopen(pymol_axis, "a+"))==NULL)
        {
            printf("\n\n** Error appending to file '%s'!",pymol_axis);
            exit(1);
        }
    }
    
   if(helix[i].residues_total>=7)
   {
       // set a flag processing new helix
       vec_flag=0;
       
      for (j=0; j<helix[i].residues_total-6; j+=3)   /* j is axis number and must be incremented by three */
      // dmf 6.28.17 exploit this j+=3 loop to handle outputing axis points for
      //     display in pymol
      {
         if((helix[i].unit_local_axis[j][0]!=-1) && (helix[i].unit_local_axis[j][1]!=-1) && (helix[i].unit_local_axis[j][2]!=-1) && (helix[i].unit_local_axis[j+3][0]!=-1) && (helix[i].unit_local_axis[j+3][1]!=-1) && (helix[i].unit_local_axis[j+3][2]!=-1))
         {
            angle=helix[i].unit_local_axis[j][0] * helix[i].unit_local_axis[j+3][0] + helix[i].unit_local_axis[j][1] * helix[i].unit_local_axis[j+3][1] + helix[i].unit_local_axis[j][2] * helix[i].unit_local_axis[j+3][2];

            helix[i].bending_angle[k]=acos(angle)*pi;

            fprintf(fpo_axis,"Bending angle between axis %d and axis %d: %f degrees\n",j,j+3,helix[i].bending_angle[k]);

            k++;   /* k is bend_number and must be incremented by one */
             
             // dmf 6.27.17 need to start assembling the PyMol script file here... (want these points as vectors)
             // - identify by helix_number + segment_number (so can create helix axis objects for each individual helix)
             // - possibly only record every third axis point (in fact, should do this, otherwise will have a ton of points)
             // - structure will be:
             //			pseudoatom pt1, pos=[x1, y1, z1]
             //			pseudoatom pt2, pos=[x2, y2, z2]
             //			distance /pt1, /pt2
             //			set dash_gap, 0
             // - use "pt1" as identifier (e.g. helix[#]pt[#]) - "i" identifies the helix, "j" identifies the segment
             // - want helix[i].origin[j][0], [1], [2] values
             // - last two lines can be added before file is closed

             // set the label for the axis point
             sprintf(pt_label,"h%dp%d",i,j);
             
             // write the coordinates of the axis point
             fprintf(fpo_pyaxis, "pseudoatom %s, pos=[ %f, %f, %f]\n", pt_label, helix[i].origin[j][0], helix[i].origin[j][1], helix[i].origin[j][2]);
             
             // write the distance statement if there's a previous point
             if (vec_flag) {
                 // may want the distance objects to have specific names (add later)
				 sprintf(dLabel,"%s_Axis",pt_label); 
                 fprintf(fpo_pyaxis, "distance %s, /%s, /%s\n", dLabel, previous_pt_label, pt_label);
                 strcpy(previous_pt_label,pt_label);
             } else {
                 strcpy(previous_pt_label,pt_label);
                 // next point, will write a distance statement
                 vec_flag = 1;
             }
         }
         else
         {
            fprintf(fpo_axis,"*** Invalid axis has been chosen ***\n");
            printf("\n\n*** Invalid axis has been chosen ***\n");
            exit(1);
         }
      }
      fprintf(fpo_axis,"\n");

      /* get maximum bending angle in helix */

      for (l=0; l<k; l++)
      {
         if(helix[i].bending_angle[l]>max_angle)
         {
            max_angle=helix[i].bending_angle[l];
         }
      }
      helix[i].max_bending_angle=max_angle;

      /* Assign geometry to helix - (K)inked or not */  

      if(max_angle>=20.0) helix[i].geometry='K';

      fprintf(fpo_axis,"Maximum bending angle: %f degrees\n\n",helix[i].max_bending_angle);
   }

   else
   {
      fprintf(fpo_axis,"Helix %d is less than 7 residues and does not have a bending angle\n\n",i);
   }

   fclose(fpo_axis);
   // dmf 7.28.17
   fclose(fpo_pyaxis); 
}    

/* ------------------------------------------------------------------------- */

/* Function to fit plane, circle, and line to the local helix origins */
void fit(int helix_number, struct HELIX *helix)
{
   /* Variables */

   int i,j,k,l;
   double rx[MAXRESIDUES], ry[MAXRESIDUES], rz[MAXRESIDUES];
   double rmp[MAXRESIDUES], rmc[MAXRESIDUES], rml[MAXRESIDUES];
   double **matp;
   double **matc;
   double **matl;       /* matp[3][3],matc[3][3],matl[2][2] */
   double rotmt[3][3];
   double **pmat;
   double **cmat;
   double **lmat;       /* pmat[3][3],cmat[3][3],lmat[2][2] */
   double v[3], ap[3], ac[3], al[2];
   double vmag;
   double x2, y2, z2, xy, xz, yz, x, y, z, x2y, xy2, y3, x3;
   double bp1,bp2,bp3,bc1,bc2,bc3,bl1,bl2;
   double sump, sumc, suml;
   double rmsdp, radcsq;
   double deno, rdeno, rnum, r;
   double pp,pl,pm,pn;
   double costheta, costheta1, theta, sintheta;
   double a1,a2,a3,b1,b2,b3;
   double rem;
   double radc, rmsdc, rmsdl, r2, ratio;
   int origins_total;
   FILE *fpo_geom;


   i=helix_number;

   origins_total=helix[i].residues_total-2;

   if(i==0)
   {
// dmf 7.25.17 - want to modify output_geom to include identifying string. 
      if((fpo_geom=fopen(output_geom, "w"))==NULL)
      {
         printf("\n\n** Error writing to file '%s'!",output_geom);
         exit(1);
      }
   }
   else   
   {
// dmf 7.25.17 - want to modify output_geom to include identifying string. 
      if((fpo_geom=fopen(output_geom, "a+"))==NULL)
      {
         printf("\n\n** Error appending to file '%s'!",output_geom);
         exit(1);
      }
   }

   fprintf(fpo_geom,"Helix Number %d\n\n",i);

   if(helix[i].residues_total>=9)
   {
      /* calloc the 3x3 and 2x2 matrices */

      matp = (double **) calloc(3,sizeof(double *));

      for(j=0; j<3; j++)
      {
         matp[j] = (double *) calloc(3,sizeof(double));
      }

      matc = (double **) calloc(3,sizeof(double *));

      for(j=0; j<3; j++)
      {
         matc[j] = (double *) calloc(3,sizeof(double));
      }

      for(j=0;j<3;j++)
      {
         for(k=0;k<3;k++)
         {
            matp[j][k]=0.0;
            matc[j][k]=0.0;
         }  
      }

      matl = (double **) calloc(2,sizeof(double *));

      for(j=0; j<2; j++)
      {
         matl[j] = (double *) calloc(2,sizeof(double));
      }

      for(j=0;j<2;j++)
      {
         for(k=0;k<2;k++)
         {
            matl[j][k]=0.0;
         }  
      }

      /* fitting least squares plane to local helix origins */

      x2=0.0;
      y2=0.0;
      z2=0.0;
      xy=0.0;
      xz=0.0;
      yz=0.0;
      x=0.0;
      y=0.0;
      z=0.0;

      for(j=0;j<origins_total;j++)
      {
         if(helix[i].origin[j][0]==-1 || helix[i].origin[j][1]==-1 || helix[i].origin[j][2]==-1)
         {
            printf("\n\n** error - invalid local helix origin used in analysis\n\n");
            exit(1);
         }

         x2 = helix[i].origin[j][0] * helix[i].origin[j][0] + x2;
         y2 = helix[i].origin[j][1] * helix[i].origin[j][1] + y2;
         z2 = helix[i].origin[j][2] * helix[i].origin[j][2] + z2;
         xy = helix[i].origin[j][0] * helix[i].origin[j][1] + xy;
         xz = helix[i].origin[j][0] * helix[i].origin[j][2] + xz;
         yz = helix[i].origin[j][1] * helix[i].origin[j][2] + yz;
         x = helix[i].origin[j][0] + x;
         y = helix[i].origin[j][1] + y;
         z = helix[i].origin[j][2] + z;
      }

      matp[0][0] = x2 + matp[0][0];
      matp[0][1] = xy + matp[0][1];
      matp[0][2] = xz + matp[0][2];
      matp[1][0] = xy + matp[1][0];
      matp[1][1] = y2 + matp[1][1];
      matp[1][2] = yz + matp[1][2];
      matp[2][0] = xz + matp[2][0];
      matp[2][1] = yz + matp[2][1];
      matp[2][2] = z2 + matp[2][2];

      pmat=matinv3(matp);

      bp1=x;
      bp2=y;
      bp3=z;

      for(j=0;j<3;j++)
      {
         ap[j]=pmat[j][0]*bp1+pmat[j][1]*bp2+pmat[j][2]*bp3;
      }

      sump=0.0;

      for(j=0;j<origins_total;j++)
      {
         rmp[j]=ap[0]*helix[i].origin[j][0] + ap[1]*helix[i].origin[j][1] + ap[2]*helix[i].origin[j][2] -1;
         sump=sump+rmp[j]*rmp[j];
      }

      if(sump<0.0) fprintf(fpo_geom,"Fit to the plane is not good\n\n");

      else
      {
         rmsdp=sqrt(sump/origins_total);
         fprintf(fpo_geom,"RMS deviation from best plane: %f\n",rmsdp);
      }

      /* converting to the normal form of the plane */
      /* deno is denominator, pp is distance of plane from XY plane */

      deno=sqrt(SQR(ap[0])+SQR(ap[1])+SQR(ap[2]));
      pp=1.0/deno;
      pl=ap[0]*pp;
      pm=ap[1]*pp;
      pn=ap[2]*pp;

      fprintf(fpo_geom,"L, M, and N of the best plane: %f, %f, %f\n",pl,pm,pn);

      /* reorientate the points so that the best fit plane coincides with the X-Y plane */
      /* l,m,n of the normal to X-Y plane are (0,0,1) */
      /* angle between normal to best fit plane & the X-Y plane */
      /* vector normal to the l,m,n of the best fit and X-Y plane... */

      for(j=0;j<3;j++)
      {
         v[j]=0.0;
      }

      v[0]=pm+v[0];
      v[1]=-pl+v[1];
      v[2]=0.0+v[2];

      vmag=sqrt(SQR(v[0])+SQR(v[1])+SQR(v[2]));
      v[0]=v[0]/vmag;
      v[1]=v[1]/vmag;
      v[2]=v[2]/vmag;

      /* the rotation matrix required to make l,m,n as 0,0,1 */

      costheta=pn;
      costheta1=1.0-pn;
      theta=acos(pn);
      sintheta=sin(theta);
      a1=v[0]*sintheta;
      a2=v[1]*sintheta;
      a3=v[2]*sintheta;
      b1=v[1]*v[2]*costheta1;
      b2=v[2]*v[0]*costheta1;
      b3=v[0]*v[1]*costheta1;

      for(j=0;j<3;j++)
      {
         for(k=0;k<3;k++)
         {
            rotmt[j][k]=0.0;
         }  
      }

      rotmt[0][0]=costheta+SQR(v[0])*costheta1+rotmt[0][0];
      rotmt[1][1]=costheta+SQR(v[1])*costheta1+rotmt[1][1];
      rotmt[2][2]=costheta+SQR(v[2])*costheta1+rotmt[2][2];
      rotmt[0][1]=b3-a3+rotmt[0][1];
      rotmt[1][0]=b3+a3+rotmt[1][0];
      rotmt[2][0]=b2-a2+rotmt[2][0];
      rotmt[0][2]=b2+a2+rotmt[0][2];
      rotmt[1][2]=b1-a1+rotmt[1][2];
      rotmt[2][1]=b1+a1+rotmt[2][1];

      /* rotating the points so that they lie in the X-Y plane */

      for(j=0;j<MAXRESIDUES;j++)
      {
         rx[j]=0.0;
         ry[j]=0.0;
         rz[j]=0.0;
      }

      for(j=0;j<origins_total;j++)
      {
         rx[j]=rotmt[0][0]*helix[i].origin[j][0]+rotmt[0][1]*helix[i].origin[j][1]+rotmt[0][2]*helix[i].origin[j][2]+rx[j];
         ry[j]=rotmt[1][0]*helix[i].origin[j][0]+rotmt[1][1]*helix[i].origin[j][1]+rotmt[1][2]*helix[i].origin[j][2]+ry[j];
         rz[j]=rotmt[2][0]*helix[i].origin[j][0]+rotmt[2][1]*helix[i].origin[j][1]+rotmt[2][2]*helix[i].origin[j][2]+rz[j];
      }

      /* translating the points so that they lie in the X-Y plane of the form ax+by=0 */

      for(j=0;j<origins_total;j++)
      {
         rx[j]=rx[j];
         ry[j]=ry[j];
         rz[j]=rz[j]-pp;
      }

      /* fitting circle to these reorientated points */
      /* equation used is (x+a)^2+(y+b)^2=r^2 */
      /* centre of circle (-a,-b) ,c=r^2-a^2-b^2 */
      /* delta=x^2+y^2+2ax+2by-c */
      /* circle finds out radius of curvature and coordinates of centre */
      /* from a set of given points using LEAST SQUARE FIT method. */

      x2 = 0.0;
      y2 = 0.0;
      x3 = 0.0;
      y3 = 0.0;
      xy2 = 0.0;
      x2y = 0.0;
      xy = 0.0;
      x = 0.0;
      y = 0.0;

      for(j=0;j<origins_total;j++)
      {
         x3 = rx[j] * rx[j] * rx[j] + x3;
         y3 = ry[j] * ry[j] * ry[j] + y3;
         x2y = rx[j] * rx[j] * ry[j] + x2y;
         xy2 = rx[j] * ry[j] * ry[j] + xy2;
         x2 = rx[j] * rx[j] + x2;
         y2 = ry[j] * ry[j] + y2;
         xy = rx[j] * ry[j] + xy;
         x = rx[j] + x;
         y = ry[j] + y;
      }

      matc[0][0]=x2+matc[0][0];
      matc[0][1]=xy+matc[0][1];
      matc[0][2]=-x+matc[0][2];
      matc[1][0]=xy+matc[1][0];
      matc[1][1]=y2+matc[1][1];
      matc[1][2]=-y+matc[1][2];
      matc[2][0]=x+matc[2][0];
      matc[2][1]=y+matc[2][1];
      matc[2][2]=-origins_total+matc[2][2];

      cmat=matinv3(matc);

      bc1=-x3-xy2;
      bc2=-x2y-y3;
      bc3=-x2-y2;

      for(j=0;j<3;j++)
      {
         ac[j]=0.0;
      }

      for(j=0;j<3;j++)
      {
         ac[j]=cmat[j][0]*bc1+cmat[j][1]*bc2+cmat[j][2]*bc3+ac[j];
      }

      ac[0]=ac[0]/2.0;
      ac[1]=ac[1]/2.0;
      ac[2]=ac[2];

      radcsq=SQR(ac[0])+SQR(ac[1])+ac[2];

      if(radcsq>0.0)
      {
         radc=sqrt(radcsq);
         fprintf(fpo_geom,"Radius of the best circle: %f\n",radc);
         fprintf(fpo_geom,"Centre of the best circle: %f, %f\n",-ac[1],-ac[2]);
      }
      else
      {
         radc=0.0;
         fprintf(fpo_geom,"Fit of the circle is not good\n\n");
      }

      sumc=0.0;

      /* dev=sqrt((x+a)^2+(y+b)^2)-r */
      /* rms dev=sqrt(dev^2/n) */

      for(j=0;j<origins_total;j++)
      {
         rem=SQR(rx[j]+ac[0])+SQR(ry[j]+ac[1]);
         rmc[j]=sqrt(rem)-radc;
         sumc=sumc+SQR(rmc[j]);
      }

      if(sumc<0.0)
      {
         fprintf(fpo_geom,"Fit of the circle is not good\n\n");
      }
      else
      {
         rmsdc=sqrt(sumc/origins_total);
         fprintf(fpo_geom,"RMS deviation from the best circle: %f\n",rmsdc);
      }

      /* Fitting line to the reorientated points - equation of line y=mx+c */

      matl[0][0]=x2+matl[0][0];
      matl[0][1]=x+matl[0][1];
      matl[1][0]=x+matl[1][0];
      matl[1][1]=origins_total+matl[1][1];

      lmat=matinv2(matl);

      bl1=xy;
      bl2=y;

      al[0]=0.0;
      al[1]=0.0;

      for(j=0;j<2;j++)
      {
         al[j]=lmat[j][0]*bl1+lmat[j][1]*bl2+al[j];
      }

      fprintf(fpo_geom,"Slope of the best line: %f\n",al[0]);
      fprintf(fpo_geom,"Intercept of the best line: %f\n",al[1]);

      suml=0.0;

      for(j=0;j<origins_total;j++)
      {
         rml[j]=ry[j]-al[0]*rx[j]-al[1];
         suml=suml+SQR(rml[j]);
      }

      if(suml<0.0) fprintf(fpo_geom,"Fit to the line is not good\n\n");

      else
      {
         rmsdl=sqrt(suml/origins_total);
         fprintf(fpo_geom,"RMS deviation from best line: %f\n",rmsdl);
      }

      /* Linear correlation coefficient */
      /* r=(n*(xy)-x*y)/sqrt((n*x2-x*x)*(n*y2-y*y)) */

      rnum=origins_total*xy-x*y;
      rdeno=(origins_total*x2-x*x)*(origins_total*y2-y*y);
      r=rnum/sqrt(rdeno);
      r2=r*r;
      fprintf(fpo_geom,"Square of linear correlation coefficient: %f\n",r2);

      /* Assign geometry to helix - (U)nknown, (C)urved or (L)inear */

      if((helix[i].max_bending_angle<20.0) && (helix[i].geometry!='K'))
      {
         ratio=rmsdl/rmsdc;

         fprintf(fpo_geom,"Ratio rmsdl/rmsdc: %f",ratio);

         if((rmsdc>1.0) && (rmsdl>1.0))
         {
            helix[i].geometry='U';
         }

         else if(ratio>1.0) helix[i].geometry='C';

         else if((ratio<=0.7) && (r2>=0.8)) helix[i].geometry='L';

         else if((ratio>0.7) && (ratio<=1.0) && (r2<=0.5)) helix[i].geometry='C';

         else if((ratio>0.7) && (ratio<=1.0) && (r2>=0.8)) helix[i].geometry='L';

         else if((ratio>0.7) && (ratio<=1.0) && (r2>0.5) && (r2<0.8)) helix[i].geometry='U';

         else helix[i].geometry='U';
      }

      fprintf(fpo_geom,"\n\n");

      destroy_matrix3(matp);
      destroy_matrix3(matc);
      destroy_matrix2(matl);
      destroy_matrix3(pmat);
      destroy_matrix3(cmat);
      destroy_matrix2(lmat);
   }

   else
   {
      fprintf(fpo_geom,"Helix %d is less than 9 residues and cannot undergo accurate line/curve fitting\n\n",i);
   }

   fclose(fpo_geom);
}

/* ------------------------------------------------------------------------- */

/* Function to calculate the inverse of a 3 x 3 matrix */
double** matinv3(double **h) 
{
   /* Variables */

   double deth, deth1, deth2, deth3;
   int i,j,k;
   double c11,c12,c13,c21,c22,c23,c31,c32,c33;
   double **s; 


   s = (double **) calloc(3,sizeof(double *));
   
   for(i=0; i<3; i++)
   {
      s[i] = (double *) calloc(3,sizeof(double));
   }

   for(j=0; j<3; j++)
   {
      for(k=0; k<3; k++)
      {
         s[j][k]=0.0;
      }
   }

   c11=h[1][1]*h[2][2]-h[2][1]*h[1][2];
   c12=h[1][0]*h[2][2]-h[2][0]*h[1][2];
   c13=h[1][0]*h[2][1]-h[2][0]*h[1][1];
   c21=h[0][1]*h[2][2]-h[0][2]*h[2][1];
   c22=h[0][0]*h[2][2]-h[0][2]*h[2][0];
   c23=h[0][0]*h[2][1]-h[2][0]*h[0][1];
   c31=h[0][1]*h[1][2]-h[0][2]*h[1][1];
   c32=h[0][0]*h[1][2]-h[0][2]*h[1][0];
   c33=h[0][0]*h[1][1]-h[0][1]*h[1][0];

   deth1=h[0][0]*c11;
   deth2=-h[0][1]*c12;
   deth3=h[0][2]*c13;
   deth=deth1+deth2+deth3;

   if(deth!=0.0)
   {
      s[0][0]=c11/deth+s[0][0];
      s[0][1]=-c21/deth+s[0][1];
      s[0][2]=c31/deth+s[0][2];
      s[1][0]=-c12/deth+s[1][0];
      s[1][1]=c22/deth+s[1][1];
      s[1][2]=-c32/deth+s[1][2];
      s[2][0]=c13/deth+s[2][0];
      s[2][1]=-c23/deth+s[2][1];
      s[2][2]=c33/deth+s[2][2];
   }
   else
   {
   printf("\n\n** ERROR ** determinant of (3x3) matrix is zero, no inverse possible!\n\n");
   }

   return s;
}

/* ------------------------------------------------------------------------- */

/* Function to calculate the inverse of a 2 x 2 matrix */
double** matinv2(double **p) 
{
   /* Variables */

   double detp, detp1, detp2;
   int i,j,k;
   double c11,c12,c21,c22;
   double **q; 


   q = (double **) calloc(2,sizeof(double *));
   
   for(i=0; i<2; i++)
   {
      q[i] = (double *) calloc(2,sizeof(double));
   }

   for(j=0; j<2; j++)
   {
      for(k=0; k<2; k++)
      {
         q[j][k]=0.0;
      }
   }

   c11=p[1][1];
   c12=p[1][0];
   c21=p[0][1];
   c22=p[0][0];

   detp1=p[0][0]*c11;
   detp2=-p[0][1]*c12;
   detp=detp1+detp2;

   if(detp!=0.0)
   {
      q[0][0]=c11/detp+q[0][0];
      q[0][1]=-c21/detp+q[0][1];
      q[1][0]=-c12/detp+q[1][0];
      q[1][1]=c22/detp+q[1][1];
   }
   else
   {
   printf("\n\n** ERROR ** determinant of (2x2) matrix is zero, no inverse possible!\n\n");
   }

   return q;
}

/* ------------------------------------------------------------------------- */

/* Free up the memory taken by the 3x3 inverse matrix (s) */
void destroy_matrix3(double **s)
{
   /* Variables */

   int i;


   for (i=0; i<3; i++)
   {
       free(s[i]);
   }
   free(s);
}

/* ------------------------------------------------------------------------- */

/* Free up the memory taken by the 2x2 inverse matrix (q) */
void destroy_matrix2(double **q)
{
   /* Variables */

   int i;


   for (i=0; i<2; i++)
   {
       free(q[i]);
   }
   free(q);
}

/* ------------------------------------------------------------------------- */

/* Function to generate a 2-D structure array of all helices against all helices in the file */
struct HELIXPAIR** neighbours(int *helices_total)
{
   /* Variables */

   int g,h,i,j,k;
   struct HELIXPAIR **helix_pair;


   /* Allocate the memory for the 2d array of HELIXPAIR dynamically
    *
    * ----no. of helices--->
    *        |
    *        |
    *  no. of helices 
    *        \/
    */

   helix_pair = (struct HELIXPAIR **) calloc(*helices_total,sizeof(struct HELIXPAIR*));
   
   for(k=0; k<*helices_total; k++)
   {
      helix_pair[k] = (struct HELIXPAIR *) calloc(*helices_total,sizeof(struct HELIXPAIR));
   }
   
   /* fill entire 2-D array with helix numbers, zeros and negative values */

   for(g=0; g<*helices_total; g++)
   {
      for(h=0; h<*helices_total; h++)
      {            
         helix_pair[g][h].helix_one=g;
         helix_pair[g][h].helix_two=h;
         helix_pair[g][h].neighbours=0;
         helix_pair[g][h].packed=0;
         helix_pair[g][h].h1_residues=-1;
         helix_pair[g][h].h2_residues=-1;
         helix_pair[g][h].h1_start=-1;
         helix_pair[g][h].h2_start=-1;
         helix_pair[g][h].h1_end=-1;
         helix_pair[g][h].h2_end=-1;
         helix_pair[g][h].vdw=0;
         helix_pair[g][h].hbond=0;
         helix_pair[g][h].electrostatic=0;
         helix_pair[g][h].covalent=0;
         helix_pair[g][h].angle1=-1.0;
         helix_pair[g][h].angle2=-1.0;
         helix_pair[g][h].distance=-1.0;
      }   
   }

return helix_pair;
}

/* ------------------------------------------------------------------------- */

/* Function to determine which helices may be packed (i.e. are neighbours) and assign that information to the HELIXPAIR structure array */
struct DISTANCE** residue_distance(int helix1, int helix2, struct HELIX *helix, struct HELIXPAIR **helix_pair)
{
   /* Variables */

   int g,h,i,j,k,m,n;
   struct DISTANCE **residue_residue;
   char temp_atom_name[5]="XXXX";
   char atom_name[5]=" CA ";


   i=helix1;
   j=helix2;

   /* Allocate the memory for the 2d array of DISTANCE dynamically
    *
    * ---- residues for helix 2 ---->
    *        |
    *        Ý
    *   residues for helix 1 
    *        \/
    */

   residue_residue = (struct DISTANCE **) calloc(helix[i].residues_total,sizeof(struct DISTANCE*));
   
   for(k=0; k<helix[i].residues_total; k++)
   {
      residue_residue[k] = (struct DISTANCE *) calloc(helix[j].residues_total,sizeof(struct DISTANCE));
   }
   
   /* fill all atom names and residue/helix numbers with junk */

   for(g=0; g<helix[i].residues_total; g++)
   {
      for(h=0; h<helix[j].residues_total; h++)
      {            
         strcpy(residue_residue[g][h].atom_name1,temp_atom_name);
         strcpy(residue_residue[g][h].atom_name2,temp_atom_name);
         residue_residue[g][h].residue_number1=-1;
         residue_residue[g][h].residue_number2=-1;
         residue_residue[g][h].helix_number1=-1;
         residue_residue[g][h].helix_number2=-1;
      }   
   }

   if(i!=j)
   {
      for(g=0; g<helix[i].residues_total; g++)
      {
         for(h=0; h<helix[j].residues_total; h++)
         {
            residue_residue[g][h].distance=sqrt(pow(helix[i].ca_coord[g][0]-helix[j].ca_coord[h][0],2.0)+pow(helix[i].ca_coord[g][1]-helix[j].ca_coord[h][1],2.0)+pow(helix[i].ca_coord[g][2]-helix[j].ca_coord[h][2],2.0));
            strcpy(residue_residue[g][h].atom_name1,atom_name);
            strcpy(residue_residue[g][h].atom_name2,atom_name);
            residue_residue[g][h].residue_number1=helix[i].residue_numbers[g];
            residue_residue[g][h].residue_number2=helix[j].residue_numbers[h];
            residue_residue[g][h].helix_number1=i;
            residue_residue[g][h].helix_number2=j;
         }
      }

      /* different helices are neighbours if they have C-alpha atoms within 10 Angstroms of each other */
      /* and the distance between their midpoint C-alphas is less than 25 Angstroms */

      m=helix[i].residues_total/2;    /* middle C-alpha of helix 1 */
      n=helix[j].residues_total/2;    /* middle C-alpha of helix 2 */

      if(residue_residue[m][n].distance<=25.0) 
      {
         for(g=0; g<helix[i].residues_total; g++)
         {
            for(h=0; h<helix[j].residues_total; h++)
            {
               if(residue_residue[g][h].distance<=10.0) 
               {
                  helix_pair[i][j].neighbours=1;     
               }
            }
         }
      }
//    printf("\n\nH%d,H%d,midpoint distance=%f Ang\n",i,j,residue_residue[m][n].distance);
   }

   return residue_residue;
}

/* ------------------------------------------------------------------------- */

/* Function that determines if neighbouring helices are packed and the nature of the interhelical contacts */
struct DISTANCE** atom_distance(int helix1, int helix2, struct HELIX *helix, struct ATOM **helix_atom, struct HELIXPAIR **helix_pair)
{
   /* Variables */

   int g,h,i,j,k,m=0,n=0,p=0,q=0,e=0,f=0;
   struct DISTANCE **atom_atom;
   float current_residue1=-1.5;
   float current_residue2=-1.5;
   float last_residue1=-1.5;
   float last_residue2=-1.5;
   float previous_residue[5];
   float first_residue1=-1.5;
   float first_residue2=-1.5;
   float atom1_vdw_rad;
   float atom1_cov_rad;
   float atom2_vdw_rad;
   float atom2_cov_rad;
   float vdw_rad_sum;
   float vdw_rad_sum2;
   float cov_rad_sum;
   float hbond_distance;
   int switch_end;
   FILE *fpo_contact;
   
   i=helix1;        
   j=helix2;

// dmf 6.27.17 
// ah! a bug - this is tested by helix number, but if helix '0' 
// comes through later in the cycle, it opens to write (see if(i==0), below) 
// rather than append. will be easy to fix if use a flag, instead.
 
// if(i==0)
// dmf 6.27.17 use a global flag
   if (new_open) 
   {
// dmf 7.25.17 - want to modify output_contact to include identifying string. 
      if((fpo_contact=fopen(output_contact, "w"))==NULL)
      {
         printf("\n\n** Error writing to file '%s'!",output_contact);
         exit(1);
      }
      new_open = 0; 
   }
   else   
   {
// dmf 7.25.17 - want to modify output_contact to include identifying string. 
      if((fpo_contact=fopen(output_contact, "a+"))==NULL)
      {
         printf("\n\n** Error appending to file '%s'!",output_contact);
         exit(1);
      }
   }

   /* Allocate the memory for the 2d array of DISTANCE dynamically
    *
    * ---- atoms for helix 2 ---->
    *        |
    *        Ý
    *   atoms for helix 1 
    *        \/
    */

   atom_atom = (struct DISTANCE **) calloc(helix[i].atoms_total,sizeof(struct DISTANCE*));
   
   for(k=0; k<helix[i].atoms_total; k++)
   {
      atom_atom[k] = (struct DISTANCE *) calloc(helix[j].atoms_total,sizeof(struct DISTANCE));
   }

   /* fill all atom names and residue/helix numbers with junk */

   for(g=0; g<helix[i].atoms_total; g++)
   {
      for(h=0; h<helix[j].atoms_total; h++)
      {            
         strcpy(atom_atom[g][h].atom_name1,"XXXX");
         strcpy(atom_atom[g][h].atom_name2,"XXXX");
         atom_atom[g][h].residue_number1=-1;
         atom_atom[g][h].residue_number2=-1;
         atom_atom[g][h].helix_number1=-1;
         atom_atom[g][h].helix_number2=-1;
      }   
   }

   for(g=0; g<5; g++)
   {
      previous_residue[g]=-1.5;
   }

   fprintf(fpo_contact,"Interhelical Contact Residues\n\n");

   for(g=0; g<helix[i].atoms_total; g++)
   {
      for(h=0; h<helix[j].atoms_total; h++)
      {
         atom_atom[g][h].distance=sqrt(pow(helix_atom[i][g].atom_coord[0]-helix_atom[j][h].atom_coord[0],2.0)+pow(helix_atom[i][g].atom_coord[1]-helix_atom[j][h].atom_coord[1],2.0)+pow(helix_atom[i][g].atom_coord[2]-helix_atom[j][h].atom_coord[2],2.0));

         strcpy(atom_atom[g][h].atom_name1,helix_atom[i][g].atom_name);

         strcpy(atom_atom[g][h].atom_name2,helix_atom[j][h].atom_name);

         /* set covalent and vdW radii */

         atom1_vdw_rad=0.0;
         atom1_cov_rad=0.0;
         atom2_vdw_rad=0.0;
         atom2_cov_rad=0.0;

         if(atom_atom[g][h].atom_name1[1]=='C')
         {
            atom1_vdw_rad=1.70;
            atom1_cov_rad=0.77;
         }
         if((atom_atom[g][h].atom_name1[1]=='O') || (atom_atom[g][h].atom_name1[0]=='O'))
         {
            atom1_vdw_rad=1.52;
            atom1_cov_rad=0.73;
         }
         if((atom_atom[g][h].atom_name1[1]=='N') || (atom_atom[g][h].atom_name1[0]=='N'))
         {
            atom1_vdw_rad=1.55;
            atom1_cov_rad=0.75;
         }
         if(atom_atom[g][h].atom_name1[1]=='S')
         {
            atom1_vdw_rad=1.80;
            atom1_cov_rad=1.02;
         }
         if(atom_atom[g][h].atom_name2[1]=='C')
         {
            atom2_vdw_rad=1.70;
            atom2_cov_rad=0.77;
         }
         if((atom_atom[g][h].atom_name2[1]=='O') || (atom_atom[g][h].atom_name2[0]=='O'))
         {
            atom2_vdw_rad=1.52;
            atom2_cov_rad=0.73;
         }
         if((atom_atom[g][h].atom_name2[1]=='N') || (atom_atom[g][h].atom_name2[0]=='N'))
         {
            atom2_vdw_rad=1.55;
            atom2_cov_rad=0.75;
         }
         if(atom_atom[g][h].atom_name2[1]=='S')
         {
            atom2_vdw_rad=1.80;
            atom2_cov_rad=1.02;
         }

         /* set hbond distances */

         hbond_distance=0.0;

         if((atom_atom[g][h].atom_name1[1]=='O' || atom_atom[g][h].atom_name1[0]=='O') && (atom_atom[g][h].atom_name2[1]=='O' || atom_atom[g][h].atom_name2[0]=='O'))
         {
            if((helix_atom[i][g].hdonor==1) || (helix_atom[j][h].hdonor==1))
            {
               hbond_distance=2.70;
            }
            else hbond_distance=0.0;
         }

         if((atom_atom[g][h].atom_name1[1]=='O' || atom_atom[g][h].atom_name1[0]=='O') && (atom_atom[g][h].atom_name2[1]=='N' || atom_atom[g][h].atom_name2[0]=='N'))
         {
            if(helix_atom[j][h].hdonor) hbond_distance=3.04; /* N is donor */

            else if((helix_atom[i][g].hdonor==1) && (helix_atom[j][h].hdonor==0)) /* O is donor and N is not */ 
            {
               hbond_distance=2.88;
            }
            else hbond_distance=0.0;
         }

         if((atom_atom[g][h].atom_name1[1]=='N' || atom_atom[g][h].atom_name1[0]=='N') && (atom_atom[g][h].atom_name2[1]=='O' || atom_atom[g][h].atom_name2[0]=='O'))
         {
            if(helix_atom[i][g].hdonor) hbond_distance=3.04; /* N is donor */

            else if((helix_atom[i][g].hdonor==0) && (helix_atom[j][h].hdonor==1)) /* O is donor and N is not */ 
            {
               hbond_distance=2.88;
            }
            else hbond_distance=0.0;
         }

         if((atom_atom[g][h].atom_name1[1]=='N' || atom_atom[g][h].atom_name1[0]=='N') && (atom_atom[g][h].atom_name2[1]=='N' || atom_atom[g][h].atom_name2[0]=='N'))
         {
            if((helix_atom[i][g].hdonor==1) || (helix_atom[j][h].hdonor==1))
            {
               hbond_distance=3.10;
            }
            else hbond_distance=0.0;
         }

         /* vdw_rad_sum used to determine residue-residue packing as according to Chothia */

         vdw_rad_sum = (atom1_vdw_rad + atom2_vdw_rad) + TOLERANCE;

         /* vdw_rad_sum2 used to detect van der Waal bonds between atoms i.e. they are within 106% of the sum of their vdW radii */
         
         vdw_rad_sum2 = (atom1_vdw_rad + atom2_vdw_rad) * TOLERANCE2;

         /* cov_rad_sum used to detect covalently bonded atoms i.e. they are within 106% of the sum of their covalent radii */

         cov_rad_sum = (atom1_cov_rad + atom2_cov_rad) * TOLERANCE2;

         atom_atom[g][h].residue_number1=helix_atom[i][g].residue_number;
         atom_atom[g][h].residue_number2=helix_atom[j][h].residue_number;

         atom_atom[g][h].helix_number1=i;
         atom_atom[g][h].helix_number2=j;

         current_residue1=atom_atom[g][h].residue_number1;
         current_residue2=atom_atom[g][h].residue_number2;

         /* residues in contact if atoms are within 0.6 A of the sum of their van der Waals' radii */

         if((atom_atom[g][h].distance <= vdw_rad_sum) && (current_residue1!=last_residue1) && (i!=j))
         {
            m++;        /* m is no. of residues of helix 1 involved in contact area */

            if(last_residue1==-1.5) first_residue1=current_residue1;   /* sequence number of first residue of helix 1 in contact area */

            last_residue1=current_residue1;                          /* sequence number of last residue of helix 1 in contact area */

#ifdef DEBUG
  	fprintf(fpo_contact,"\t\tTesting g %d and h %d\n", g, h);
#endif

            fprintf(fpo_contact,"Helix %d residue: %.2f\n",i,last_residue1);
         }

         if((atom_atom[g][h].distance <= vdw_rad_sum) && (current_residue2!=previous_residue[0]) && (current_residue2!=previous_residue[1]) && (current_residue2!=previous_residue[2]) && (current_residue2!=previous_residue[3]) && (current_residue2!=previous_residue[4]) && (i!=j))
         {
            n++;        /* n is no. of residues of helix 2 involved in contact area */

            if(previous_residue[0]==-1.5)           /* happens on first occasion only */
            {
               first_residue2=current_residue2;   /* sequence number of first residue of helix 2 in contact area */

               last_residue2=current_residue2;    /* sequence number of last residue of helix 2 in contact area */
            }

            if(current_residue2<first_residue2) first_residue2=current_residue2;  /* make first_residue2 the lowest residue number */

            if(current_residue2>last_residue2) last_residue2=current_residue2;    /* make last_residue2 the highest residue number */

            for(k=4; k>0; k--)         /* update backcatalogue of past 5 residues of helix 2 to prevent recording same residue twice */
            {
               previous_residue[k]=previous_residue[k-1];
            }

            previous_residue[0]=current_residue2;

#ifdef DEBUG
  	fprintf(fpo_contact,"\t\tTesting g %d and h %d\n", g, h);
#endif
            
            fprintf(fpo_contact,"Helix %d residue: %.2f\n",j,previous_residue[0]);
         }

         if((atom_atom[g][h].distance <= cov_rad_sum) && (i!=j))
         {
            q++;        /* q is no. of covalent contacts between two helices */ 
         }

         else if((atom_atom[g][h].distance <= 5.0) && (helix_atom[i][g].charge+helix_atom[j][h].charge==3) && (i!=j))
         {
            e++;        /* e is no. of electrostatic interhelical contacts */
         }

         /* Deemed to be H-bonded if within 106% of the appropriate H-bond distance */

         else if((atom_atom[g][h].distance <= hbond_distance * TOLERANCE2) && (i!=j))
         {
            f++;        /* f is no. of H-bonds between the two helices */
         }

         else if((atom_atom[g][h].distance <= vdw_rad_sum2) && (i!=j))
         {
            p++;        /* p is no. of vdW contacts between two helices */
         }
      }
   }

   /* two helices are packed if each has at least three residues forming contacts to the other */

   if((m>=MIN_CONTACT_RESIDUES) && (n>=MIN_CONTACT_RESIDUES)) 
   {
      helix_pair[i][j].packed=1;
   }

   fprintf(fpo_contact,"Helix %d first residue: %.2f\t",i,first_residue1);
   fprintf(fpo_contact,"Helix %d last residue: %.2f\n",i,last_residue1);
   fprintf(fpo_contact,"Helix %d first residue: %.2f\t",j,first_residue2);
   fprintf(fpo_contact,"Helix %d last residue: %.2f\n\n",j,last_residue2);

   /* record the number of contacting residues that are involved in the packed pair */

   helix_pair[i][j].h1_residues=m;
   helix_pair[i][j].h2_residues=n;

   /* find the start and end points of the contact zone for each helix */

   for(g=0; g<helix[i].residues_total; g++)
   {
      if(helix[i].residue_numbers[g]==first_residue1) helix_pair[i][j].h1_start=g;      
      if(helix[i].residue_numbers[g]==last_residue1) helix_pair[i][j].h1_end=g;      
   }

   for(g=0; g<helix[j].residues_total; g++)
   {
      if(helix[j].residue_numbers[g]==first_residue2) helix_pair[i][j].h2_start=g;      
      if(helix[j].residue_numbers[g]==last_residue2) helix_pair[i][j].h2_end=g;      
   }

   /* if the end of the contact zone is before the start (which should never happen) */
   /* but does if the helix is numbered from high to low (i.e. the wrong way round)  */
   /* then switch the end point number with the start point number and vice versa    */

   if(helix_pair[i][j].h1_end<helix_pair[i][j].h1_start)
   {
      switch_end=helix_pair[i][j].h1_end;

      helix_pair[i][j].h1_end=helix_pair[i][j].h1_start;

      helix_pair[i][j].h1_start=switch_end;
   }

   if(helix_pair[i][j].h2_end<helix_pair[i][j].h2_start)
   {
      switch_end=helix_pair[i][j].h2_end;

      helix_pair[i][j].h2_end=helix_pair[i][j].h2_start;

      helix_pair[i][j].h2_start=switch_end;
   }

   /* record interhelical bond types for that pair */

   helix_pair[i][j].vdw=p;
   helix_pair[i][j].covalent=q;
   helix_pair[i][j].electrostatic=e;
   helix_pair[i][j].hbond=f;

   fclose(fpo_contact);

   return atom_atom;
}

/* ------------------------------------------------------------------------- */

/* Free up the memory taken by residue_residue structure */
void destroy_residue_residue(struct DISTANCE **residue_residue, struct HELIX *helix, int helix_number)
{
   /* Deallocate the memory for the 2d array dynamically
    *
    * ----residues for helix 2--->
    *        |
    *  residues for helix 1
    *        \/
    */

   int i;
   int j;


   j=helix_number;

   for (i=0; i<helix[j].residues_total; i++)
   {
      free(residue_residue[i]);
   }
   free(residue_residue);
}

/* ------------------------------------------------------------------------- */

/* Free up the memory taken by atom_atom structure */
void destroy_atom_atom(struct DISTANCE **atom_atom, struct HELIX *helix, int helix_number)
{
   /* Deallocate the memory for the 2d array dynamically
    *
    * ---atoms for helix 2 --->
    *        |
    *  atoms for helix 1
    *        \/
    */

   int i;
   int j;


   j=helix_number;

   for (i=0; i<helix[j].atoms_total; i++)
   {
      free(atom_atom[i]);
   }
   free(atom_atom);
}

/* ------------------------------------------------------------------------- */

/* Function to fill up two point structures and two vector structures with information from two packed helices */
void two_helix_all_vectors(int helix_A, int helix_B, struct HELIX *helix, struct HELIXPAIR **helix_pair)
{
   /* Variables */

   int g,h,i,j;
   static POINT A;     /* start point of helix A vector                */
   static POINT B;     /* start point of helix B vector                */
   static VECTOR dA;   /* vector of helix A                            */
   static VECTOR dB;   /* vector of helix B                            */
   static POINT pA;    /* point of closest approach on line of helix A */
   static POINT pB;    /* point of closest approach on line of helix B */
   int rval;           /* 0 if parallel; 1 if lines intersect; 2 if lines are skew (as expected) */  
   double totalx,totaly,totalz;
   double averagex,averagey,averagez;
   double contact_vector[3];          /* vector of the closest approach from Helix B to Helix A */
   double normalA[3];                 
   double normalB[3];                
   double dot_normal;
   double mag_normalA;
   double mag_normalB;
   double cos_theta;
   double theta;
   double pi;
   double angle;
   double volume;


   pi=180.0/acos(-1.0);

   i=helix_A;
   j=helix_B;

   if((helix[i].residues_total<4) || (helix[j].residues_total<4))
   {
      printf("\n\n** Error ** Packed helix is less than 4 residues!\n\n");
      exit(1);
   }
   else
   {
      /* assign two helix origins from the two helices to the two start point structures A and B */

      h=(helix[i].residues_total-3)/2;   /* h = halfway along helix A */

      A.px=helix[i].origin[h][0];
      A.py=helix[i].origin[h][1];
      A.pz=helix[i].origin[h][2];

      h=(helix[j].residues_total-3)/2;   /* h = halfway along helix B */

      B.px=helix[j].origin[h][0];
      B.py=helix[j].origin[h][1];
      B.pz=helix[j].origin[h][2];

      /* average all the helix A axis vectors */

      totalx=0.0;
      totaly=0.0;
      totalz=0.0;

      for(g=0;g<helix[i].residues_total-3; g++)
      {
         totalx=helix[i].unit_local_axis[g][0]+totalx;
         totaly=helix[i].unit_local_axis[g][1]+totaly;
         totalz=helix[i].unit_local_axis[g][2]+totalz;
      }

      averagex=totalx/(helix[i].residues_total-3);
      averagey=totaly/(helix[i].residues_total-3);
      averagez=totalz/(helix[i].residues_total-3);

      /* assign the average vector of helix A to the vector structure dA */

      dA.dx=averagex;
      dA.dy=averagey;
      dA.dz=averagez;

      /* average all the helix B axis vectors */

      totalx=0.0;
      totaly=0.0;
      totalz=0.0;

      for(g=0;g<helix[j].residues_total-3; g++)
      {
         totalx=helix[j].unit_local_axis[g][0]+totalx;
         totaly=helix[j].unit_local_axis[g][1]+totaly;
         totalz=helix[j].unit_local_axis[g][2]+totalz;
      }

      averagex=totalx/(helix[j].residues_total-3);
      averagey=totaly/(helix[j].residues_total-3);
      averagez=totalz/(helix[j].residues_total-3);

      /* assign the average vector of helix B to the vector structure dB */

      dB.dx=averagex;
      dB.dy=averagey;
      dB.dz=averagez;

      // dmf 6.28.17 so A, B have an origin, dA,dB have an average vector, and
      // function returns pA, pB which are points of closest approach.
      rval=line_line_closest_points3d(&pA, &pB, &A, &dA, &B, &dB);

      /* smallest distance between the two helix axes i.e. length of line of closest approach */
      // dmf 6.28.17 Determine this in the 'contact region' routine, below.
      // helix_pair[i][j].distance=sqrt(pow(pB.px-pA.px,2.0)+pow(pB.py-pA.py,2.0)+pow(pB.pz-pA.pz,2.0));
       
      /* the vector of the closest approach from helix axis B to helix axis A */

      contact_vector[0]=pB.px-pA.px;
      contact_vector[1]=pB.py-pA.py;
      contact_vector[2]=pB.pz-pA.pz;

      /* normal vectors are cross-product of average axis vector and contact vector */

      normalA[0]=(dA.dy * contact_vector[2]) - (dA.dz * contact_vector[1]);
      normalA[1]=(dA.dz * contact_vector[0]) - (dA.dx * contact_vector[2]);
      normalA[2]=(dA.dx * contact_vector[1]) - (dA.dy * contact_vector[0]);

      normalB[0]=(dB.dy * contact_vector[2]) - (dB.dz * contact_vector[1]);
      normalB[1]=(dB.dz * contact_vector[0]) - (dB.dx * contact_vector[2]);
      normalB[2]=(dB.dx * contact_vector[1]) - (dB.dy * contact_vector[0]);

      /* get dot-product of the two normal vectors */

      dot_normal = normalA[0]*normalB[0] + normalA[1]*normalB[1] + normalA[2]*normalB[2];

      /*     nA.nB = ÝnAÝ * ÝnBÝ * cos_theta     */

      mag_normalA=sqrt(SQR(normalA[0])+SQR(normalA[1])+SQR(normalA[2]));

      mag_normalB=sqrt(SQR(normalB[0])+SQR(normalB[1])+SQR(normalB[2]));

      cos_theta=dot_normal/(mag_normalA * mag_normalB);

      angle=acos(cos_theta)*pi;

      /* V = (unit_helix_B_vector^contact_vector).unit_helix_A_vector */
      /* this is the same as dot product of normalB and unit_helix_A_vector */

      volume=(normalB[0]*dA.dx + normalB[1]*dA.dy + normalB[2]*dA.dz);

      if(volume<0) angle=-angle;
       
      // dmf 6.29.17 This is really all that this routine is trying to determine:
      helix_pair[i][j].angle1=angle;
   }
}

/* ------------------------------------------------------------------------- */

/* Function to fill up two point structures and two vector structures with information from two packed helices */
void two_helix_contact_vectors(int helix_A, int helix_B, struct HELIX *helix, struct HELIXPAIR **helix_pair)
{
   /* Variables */

   int g,h,i,j,k,m,n;
   static POINT A;     /* start point of helix A vector                */
   static POINT B;     /* start point of helix B vector                */
   static VECTOR dA;   /* vector of helix A                            */
   static VECTOR dB;   /* vector of helix B                            */
   static POINT pA;    /* point of closest approach on line of helix A */
   static POINT pB;    /* point of closest approach on line of helix B */
   int rval;           /* 0 if parallel; 1 if lines intersect; 2 if lines are skew (as expected) */  
   int hA_start;       /* first axis in contact area of helix A        */
   int hA_end;         /* last axis in contact area of helix A         */
   int hB_start;       /* first axis in contact area of helix B        */
   int hB_end;         /* last axis in contact area of helix B         */
   double totalx,totaly,totalz;
   double averagex,averagey,averagez;
   double contact_vector[3];          /* vector of the closest approach from Helix B to Helix A */
   double normalA[3];                 
   double normalB[3];                
   double dot_normal;
   double mag_normalA;
   double mag_normalB;
   double cos_theta;
   double theta;
   double pi;
   double angle;
   double volume;
    
   // dmf 6.28.17
    char cA_label[20], cB_label[20], cD_label[20]; 
    FILE *fpo_pyaxis;
    static LINESEGMENT Alimits;    // helix A axial start and end points
    static LINESEGMENT Blimits;    // helix B axial start and end points
    float rdist;  // return value for linesegment distance routine

// dmf 7.25.17 - want to modify pymol_axis to include identifying string. 
    if((fpo_pyaxis=fopen(pymol_axis, "a+"))==NULL)
    {
        printf("\n\n** Error appending to file '%s'!",pymol_axis);
        exit(1);
    }
    
   pi=180.0/acos(-1.0);

   i=helix_A;
   j=helix_B;

   if((helix[i].residues_total<4) || (helix[j].residues_total<4))
   {
      printf("\n\n** Error ** Packed helix is less than 4 residues!\n\n");
      exit(1);
   }
   else
   {
      // dmf 6.29.17 these are the axis segment identifiers at the start and
      // end of the helix contacts (i.e. they are integers)
      hA_start=helix_pair[i][j].h1_start;
      hA_end=helix_pair[i][j].h1_end;

      hB_start=helix_pair[i][j].h2_start;
      hB_end=helix_pair[i][j].h2_end;

      /* assign the two middle helix origins from the two helices to the two start point structures A and B */

      h=(hA_start+hA_end)/2;   /* h = halfway along helix A contact region */

      if(h>helix[i].residues_total-3)
      {
         h=helix[i].residues_total-3;  /* if h is bigger than last origin, make it equal last origin */
      }

      A.px=helix[i].origin[h][0];
      A.py=helix[i].origin[h][1];
      A.pz=helix[i].origin[h][2];

      k=(hB_start+hB_end)/2;   /* k = halfway along helix B contact region */

      if(k>helix[j].residues_total-3)
      {
         k=helix[j].residues_total-3;  /* if k is bigger than last origin, make it equal last origin */
      }

      B.px=helix[j].origin[k][0];
      B.py=helix[j].origin[k][1];
      B.pz=helix[j].origin[k][2];

       // dmf 6.29.17 changed criteria to <> 30 residues
      if(hA_end-hA_start+1<33)   /* if contact region is less than 31 residues... */
      {
         /* get the helix A axis vectors from the contact zone and average */

         if(hA_start>(helix[i].residues_total-4))  /* if first axis is bigger than final axis make it equal final axis */
         {
            hA_start=helix[i].residues_total-4;
         }

         if(hA_end>(helix[i].residues_total-4))    /* if last axis is bigger than final axis make it equal final axis */
         {
            hA_end=helix[i].residues_total-4;
         }

         totalx=0.0;
         totaly=0.0;
         totalz=0.0;

         for(g=hA_start;g<=hA_end;g++)
         { 
            totalx=helix[i].unit_local_axis[g][0]+totalx;
            totaly=helix[i].unit_local_axis[g][1]+totaly;
            totalz=helix[i].unit_local_axis[g][2]+totalz;
         }

         averagex=totalx/(hA_end-hA_start+1);
         averagey=totaly/(hA_end-hA_start+1);
         averagez=totalz/(hA_end-hA_start+1);

         /* assign the average vector of first helix to the vector structure dA */

         dA.dx=averagex;
         dA.dy=averagey;
         dA.dz=averagez;
          
          // dmf 6.29.17 Have an origin (at the center) and have an axis vector; should
          // now be able to define a 'start' and 'end' point as origin +/- vector*((hA_end-hAstart+1)/2)
          // which will yield points on the averaged helix axis vector. This can be used to
          // determine the segment-wise point of close approach (which can be different from the perpendicular)
          // Alternatively, just grab the coordinates at hA_start and hA_end to define a helix vector.
          Alimits.P0.px=A.px-(totalx/2);
          Alimits.P0.py=A.py-(totaly/2);
          Alimits.P0.pz=A.pz-(totalz/2);
          Alimits.P1.px=A.px+(totalx/2);
          Alimits.P1.py=A.py+(totaly/2);
          Alimits.P1.pz=A.pz+(totalz/2);
#ifdef DEBUG
          printf("\n testing start and end points helix %d\n",i);
          printf("Center: A %f %f %f\n",A.px, A.py, A.pz);
          printf("Vector: dA %f %f %f\n", dA.dx, dA.dy, dA.dz);
          printf("Start: Ai %f %f %f\n", Alimits.P0.px, Alimits.P0.py, Alimits.P0.pz);
          printf("Endpoint: Af %f %f %f\n", Alimits.P1.px, Alimits.P1.py, Alimits.P1.pz);
#endif

      }

       // dmf 6.29.17 changed criteria to <> 30 residues
      else if(hA_end-hA_start+1>=33)  /* if the contact zone is bigger than 30 residues... (to be used mainly for coiled coils) */
      {
         totalx=0.0;
         totaly=0.0;
         totalz=0.0;

         for(m=h-12;m<h+13;m++)   /* get middle 25 axes only */
         {
            if(m<0)
            {
               printf("\n\n** Error in contact axes - less than zero - helix %d **\n",i);
               exit(1);
            }

            totalx=helix[i].unit_local_axis[m][0]+totalx;
            totaly=helix[i].unit_local_axis[m][1]+totaly;
            totalz=helix[i].unit_local_axis[m][2]+totalz;

            if((helix[i].unit_local_axis[m][0]==-1)||(helix[i].unit_local_axis[m][1]==-1)||(helix[i].unit_local_axis[m][2]==-1))
            {
               printf("\n\n** Error in contact axes - over limit - helix %d **\n",i);
               exit(1);
            }
         }

         averagex=totalx/25;
         averagey=totaly/25;
         averagez=totalz/25;

         /* assign the average vector of first helix to the vector structure dA */

         dA.dx=averagex;
         dA.dy=averagey;
         dA.dz=averagez;
      }

       // dmf 6.29.17 changed criteria to <> 30 residues
      if(hB_end-hB_start+1<33)
      {
         /* get the helix B axis vectors from the contact zone and average */

         if(hB_start>(helix[j].residues_total-4))
         {
            hB_start=helix[j].residues_total-4;
         }

         if(hB_end>(helix[j].residues_total-4))
         {
            hB_end=helix[j].residues_total-4;
         }

         totalx=0.0;
         totaly=0.0;
         totalz=0.0;

         for(g=hB_start;g<=hB_end;g++)
         {
            totalx=helix[j].unit_local_axis[g][0]+totalx;
            totaly=helix[j].unit_local_axis[g][1]+totaly;
            totalz=helix[j].unit_local_axis[g][2]+totalz;
         }

         averagex=totalx/(hB_end-hB_start+1);
         averagey=totaly/(hB_end-hB_start+1);
         averagez=totalz/(hB_end-hB_start+1);

         /* assign the average vector of second helix to the vector structure dB */

         dB.dx=averagex;
         dB.dy=averagey;
         dB.dz=averagez;
          
          // dmf 6.29.17 Have an origin (at the center) and have an axis vector; should
          // now be able to define a 'start' and 'end' point as origin +/- vector*((hA_end-hAstart+1)/2)
          // which will yield points on the averaged helix axis vector. This can be used to
          // determine the segment-wise point of close approach (which can be different from the perpendicular)
          // Alternatively, just grab the coordinates at hA_start and hA_end to define a helix vector.
          Blimits.P0.px=B.px-(totalx/2);
          Blimits.P0.py=B.py-(totaly/2);
          Blimits.P0.pz=B.pz-(totalz/2);
          Blimits.P1.px=B.px+(totalx/2);
          Blimits.P1.py=B.py+(totaly/2);
          Blimits.P1.pz=B.pz+(totalz/2);
#ifdef DEBUG
          printf("\n testing start and end points helix %d\n",j);
          printf("Center: B %f %f %f\n",B.px, B.py, B.pz);
          printf("Vector: dB %f %f %f\n", dB.dx, dB.dy, dB.dz);
          printf("Start: Bi %f %f %f\n", Blimits.P0.px, Blimits.P0.py, Blimits.P0.pz);
          printf("Endpoint: Bf %f %f %f\n", Blimits.P1.px, Blimits.P1.py, Blimits.P1.pz);
#endif

      }

       // dmf 6.29.17 changed criteria to <> 30 residues
      else if(hB_end-hB_start+1>=33)   /* for coiled coils - if contact zone is equal to 30 residues or more */
      {
         totalx=0.0;
         totaly=0.0;
         totalz=0.0;

         for(n=k-12;n<k+13;n++)   /* get middle 25 axes only */
         {
            if(n<0)
            {
               printf("\n\n** Error in contact axes - less than zero - helix %d **\n",j);
               exit(1);
            }

            totalx=helix[j].unit_local_axis[n][0]+totalx;
            totaly=helix[j].unit_local_axis[n][1]+totaly;
            totalz=helix[j].unit_local_axis[n][2]+totalz;

            if((helix[j].unit_local_axis[n][0]==-1)||(helix[j].unit_local_axis[n][1]==-1)||(helix[j].unit_local_axis[n][2]==-1))
            {
               printf("\n\n** Error in contact axes - over limit - helix %d **\n",j);
               exit(1);
            }
         }

         averagex=totalx/25;
         averagey=totaly/25;
         averagez=totalz/25;

         /* assign the average vector of second helix to the vector structure dB */

         dB.dx=averagex;
         dB.dy=averagey;
         dB.dz=averagez;
      }

      // dmf 6.28.17 so A, B have an origin, dA,dB have an average vector, and
      // function returns pA, pB which are points of closest approach
      rval=line_line_closest_points3d(&pA, &pB, &A, &dA, &B, &dB);
       
       // determine the crossing angle...
       
       /* the vector of the closest approach from helix axis B to helix axis A */
       
       contact_vector[0]=pB.px-pA.px;
       contact_vector[1]=pB.py-pA.py;
       contact_vector[2]=pB.pz-pA.pz;
       
       /* normal vectors are cross-product of average axis vector and contact vector */
       
       normalA[0]=(dA.dy * contact_vector[2]) - (dA.dz * contact_vector[1]);
       normalA[1]=(dA.dz * contact_vector[0]) - (dA.dx * contact_vector[2]);
       normalA[2]=(dA.dx * contact_vector[1]) - (dA.dy * contact_vector[0]);
       
       normalB[0]=(dB.dy * contact_vector[2]) - (dB.dz * contact_vector[1]);
       normalB[1]=(dB.dz * contact_vector[0]) - (dB.dx * contact_vector[2]);
       normalB[2]=(dB.dx * contact_vector[1]) - (dB.dy * contact_vector[0]);
       
       /* get dot-product of the two normal vectors */
       
       dot_normal = normalA[0]*normalB[0] + normalA[1]*normalB[1] + normalA[2]*normalB[2];
       
       /*     nA.nB = ÝnAÝ * ÝnBÝ * cos_theta     */
       
       mag_normalA=sqrt(SQR(normalA[0])+SQR(normalA[1])+SQR(normalA[2]));
       
       mag_normalB=sqrt(SQR(normalB[0])+SQR(normalB[1])+SQR(normalB[2]));
       
       cos_theta=dot_normal/(mag_normalA * mag_normalB);
       
       angle=acos(cos_theta)*pi;
       
       /* V = (unit_helix_B_vector^contact_vector).unit_helix_A_vector */
       /* this is the same as dot product of normalB and unit_helix_A_vector */
       
       volume=(normalB[0]*dA.dx + normalB[1]*dA.dy + normalB[2]*dA.dz);
       
       if(volume<0) angle=-angle;
       
       // dmf 6.29.17 This is really all this routine is trying to determine:
       helix_pair[i][j].angle2=angle;

       // now deal with the distance...
       
       // dmf 6.29.17 Important to realize that this value is not necessarily physical.
       // The computation of the helix crossing angle requires the approach taken with the
       // line_line_closest_points3d() call, but the resulting 'closest' point may be beyond
       // the limits of the helix.
       
       /* smallest distance between the two helix axes i.e. length of line of closest approach */
       // dmf 7.12.17 recalcualate this below (keep here for debug comparison below, als)
       helix_pair[i][j].distance=sqrt(pow(pB.px-pA.px,2.0)+pow(pB.py-pA.py,2.0)+pow(pB.pz-pA.pz,2.0));
       
       // BUG2 - GETTING contact vectors well displaced from position of helices... 
       // dmf 6.28.17 also, should output the points pA and pB for display in pymol
#ifdef DEBUG
	printf("\n about to write a contact pair vector");
	printf("\n for helices %d and %d",i, j); 
	printf("\n here are A %f %f %f", A.px, A.py, A.pz);
	printf("\n and B %f %f %f", B.px, B.py, B.pz); 
        printf("\n and here are the positions of the 'contact' vector \n");
	printf("%f %f %f\n",pA.px, pA.py, pA.pz);
	printf("%f %f %f\n",pB.px, pB.py, pB.pz); 
	printf("  \n");
#endif
       
       // check and redetermine the distance within the domain of the physical helices
       // using the values of the inital and final points Ai,Af & Bi,Bf along the helix axes.
       
      // function returns pA, pB which are points of closest approach within the segment
       rdist=segment_segment_closest_points3d(&pA, &pB, &Alimits, &Blimits);
       
#ifdef DEBUG
       printf("coming out of seg_seg routine \n");
       printf("rdist is %f",rdist);
       printf(" compare with previously calculated distance %f\n",helix_pair[i][j].distance);
       printf("segment point A %f, %f, %f\n",pA.px, pA.py, pA.pz);
       printf("segment point B %f, %f, %f\n",pB.px, pB.py, pB.pz);
#endif
       
       sprintf(cA_label,"%dto%d",i,j);
       sprintf(cB_label,"%dto%d",j,i);
	   sprintf(cD_label,"Contact_%dto%d",i,j);
       fprintf(fpo_pyaxis,"pseudoatom %s, pos=[%f, %f, %f]\n",cA_label,pA.px,pA.py,pA.pz);
       fprintf(fpo_pyaxis,"pseudoatom %s, pos=[%f, %f, %f]\n",cB_label,pB.px,pB.py,pB.pz);
       fprintf(fpo_pyaxis,"distance %s, /%s, /%s\n",cD_label, cA_label,cB_label);
       
       /* smallest distance between the two helix axes i.e. length of line of closest approach */
       // dmf 7.12.17 added this so that helix_packing_pair.txt output is consistent
       helix_pair[i][j].distance=sqrt(pow(pB.px-pA.px,2.0)+pow(pB.py-pA.py,2.0)+pow(pB.pz-pA.pz,2.0));
       
   }
   // dmf 7.28.17
   fclose(fpo_pyaxis); 
}

/* ------------------------------------------------------------------------- */

/* Free up the memory taken by helix_atom structure */
void destroy_helix_atom(struct ATOM **helix_atom, int *helices_total)
{
   /* Deallocate the memory for the 2d array dynamically
    *
    * ---max atoms per helix (1000-2000)-->
    *        |
    *  number of helices
    *        \/
    */

   int i;


   for (i=0; i<*helices_total; i++)
   {
      free(helix_atom[i]);
   }
   free(helix_atom);
}

/* ------------------------------------------------------------------------- */

/* Free up the memory taken by helix_pair structure */
void destroy_helix_pair(struct HELIXPAIR **helix_pair, int *helices_total)
{
   /* Deallocate the memory for the 2d array dynamically
    *
    * ---number of helices--->
    *        |
    *  number of helices
    *        \/
    */

   int i;


   for (i=0; i<*helices_total; i++)
   {
      free(helix_pair[i]);
   }
   free(helix_pair);
}

// dmf 7.29.17
void create_filenames(char *pdb_id) {
    /* create the output filenames with pdb_identifier appended */
    
    /* variables declared globally */
    strcpy(output_shape, pdb_id);
    strcpy(output_packing, pdb_id);
    strcpy(output_helices, pdb_id);
    strcpy(pymol_axis, pdb_id);
    strcpy(output_axis, pdb_id);
    strcpy(output_contact, pdb_id);
    strcpy(output_geom, pdb_id);
    /* */
    strcat(output_shape, "_helix_shape.txt");
    strcat(output_packing, "_helix_packing_pair.txt");
    strcat(output_helices, "_helices.txt");
    strcat(pymol_axis, "_axis.py");
    strcat(output_axis, "_axis.txt");
    strcat(output_contact, "_contact.txt");
    strcat(output_geom, "_geom.txt");
    
    printf("Helix Output Files: %s, %s, %s, %s\n", output_helices, output_axis, output_geom, output_contact);
    printf("SSE Appended Output Files: %s, %s\n",output_packing,output_shape);

}

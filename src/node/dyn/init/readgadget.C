
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Convert GADGET unformatted snapshot output:
////
////              (N, ndim, time,
////               mass[i], i = 1,...,N,
////               pos[i],  i = 1,...,N,
////               vel[i],  i = 1,...,N)
////
//// into a Starlab snapshot.
////
//// Usage:  readgadget [OPTIONS]
////
//// Options:
////              -i    number the particles sequentially [don't number]
////              -w    Write GADGET snapshot insted of reading it
////
//// Written by Piero Spinnato and Simon Portegies Zwart.
////
//// Report bugs to starlab@sns.ias.edu.

//   First version: Piero Spinnato, April 2002
//   converted to starlab: SPZ in Amsterdam May 2002

#include "dyn.h"

#include "GDT_allvars.h"
#include "GDT_proto.h"
#include "GDT_allvars1.h"

#ifdef TOOLBOX

local dyn* read_gadget(bool i_flag) {

#define SKIP my_fread(&blklen,sizeof(int4byte),1,fd);
  FILE *fd;
  int   i,k,massflag,count;
  float dummy[3];
  vec vdummy;
  int   pc,type ;
  int4byte  intdummy, blklen;
  double u_init;

  int N_skip;
  dyn *root, *by, *bo;

  // Create root node.

  root = new dyn();
  if (i_flag) root->set_label(0);


  fd=stdin;
    
  SKIP; 
  if(blklen!=256)
    {
      fprintf(stderr, "incorrect header format (1)\n");
      endrun(888);
    }
  my_fread(&header1,sizeof(header1),1,fd);
  SKIP;
  if(blklen!=256)
    {
      fprintf(stderr, "incorrect header format (2)\n");
      endrun(889);
    }

  All.TotN_gas  = N_gas  = header1.npart[0];
  All.TotN_halo = header1.npart[1];
  All.TotN_disk = header1.npart[2];
  All.TotN_bulge= header1.npart[3];
  All.TotN_stars= header1.npart[4];
          
  All.Time = All.TimeBegin = header1.time;

    
  int n; n = All.TotN_stars; PRL(n);        
  int ndim; ndim = 3; PRL(ndim);            
  real time; time = All.Time; PRL(time);   

  for(i=0, massflag=0;i<5;i++)
    {
      All.MassTable[i]= header1.mass[i];
      if(All.MassTable[i]==0 && header1.npart[i]>0) {massflag=1;}
    }
  fprintf(stderr, "\nN_sph: %d\nN_halo: %d\nN_disk: %d\nN_bulge: 
              %d\nN_stars: %d\n", All.TotN_gas, All.TotN_halo, All.TotN_disk, 
	  All.TotN_bulge, All.TotN_stars);
       
  N_skip = All.TotN_gas + All.TotN_halo + All.TotN_disk + All.TotN_bulge;

  root->set_system_time(time);
    
  // Create first daughter node.

  bo = new dyn();
  root->set_oldest_daughter(bo);
  bo->set_parent(root);
  if (i_flag) bo->set_label(1);

  // Create other daughter nodes.

  for (int i = 1; i < n; i++) {
    by = new dyn();
    if (i_flag) by->set_label(i+1);
    by->set_parent(root);
    bo->set_younger_sister(by);
    by->set_elder_sister(bo);
    by->set_mass(bo->get_mass());
    bo = by;
  }


  SKIP;
  for (i = 1; i <= N_skip; i++) { my_fread(&dummy[0],sizeof(float),3,fd); }
      
  for_all_daughters(dyn, root, b) {

    my_fread(&dummy[0],sizeof(float),3,fd);
    vdummy[0] = dummy[0]; vdummy[1] = dummy[1]; vdummy[2] = dummy[2];
    b->set_pos(vdummy);
  }
  SKIP;


  SKIP;
  for (i = 1; i <= N_skip; i++) { my_fread(&dummy[0],sizeof(float),3,fd); }

  for_all_daughters(dyn, root, b){

    my_fread(&dummy[0],sizeof(float),3,fd);
    vdummy[0] = dummy[0]; vdummy[1] = dummy[1]; vdummy[2] = dummy[2];
    b->set_vel(vdummy);
  }
  SKIP;


  SKIP;
  for (i = 1; i <= N_skip; i++) { my_fread(&intdummy, sizeof(int4byte), 1, fd); }
  for_all_daughters(dyn, root, b)
    {
      my_fread(&intdummy, sizeof(int4byte), 1, fd);
    }
  SKIP;

    
  real total_mass = 0;
    
  if(massflag) { SKIP; }

  for(type=0, count=1; type<4; type++)
    {
      if(All.MassTable[type]==0 && header1.npart[type]>0)
	{
	  for(i=1;i<=header1.npart[type];i++)
	    {
	      my_fread(&dummy[0],sizeof(float),1,fd);
	    } } }

  for_all_daughters(dyn, root, b) {
            
    if(All.MassTable[4]==0 && header1.npart[4]>0)
      {
	my_fread(&dummy[0],sizeof(float),1,fd);
      }
    else { dummy[0] = All.MassTable[4]; }

    b->set_mass(dummy[0]);
    total_mass += dummy[0];
  }

  root->set_mass(total_mass);

  return root;
}

local void write_gadget() {

  int   i,k,massflag,count, idum;
  vec dummy, du0;
  int   pc, type ;
  int4byte  intdummy, blklen;
  double u_init;


  dyn *b;
  b = get_dyn(cin);
  b->flatten_node();
  
  // for safetly check number of leaves.
  int n=0;
  for_all_daughters(dyn, b, bi) {
    n++;
  }
  
  //start writing GADGET input file (code adapted from S-Gadget routines)

  All.PartAllocFactor = 1;  /* needed for memory allocation */ 
    
  All.TotN_gas  = N_gas  = 0;
  All.TotN_halo = 0;
  All.TotN_disk = 0;
  All.TotN_bulge= 0;
  All.TotN_stars= n;
      
  All.Time = All.TimeBegin = b->get_system_time();

  All.MassTable[0]= 1;
  All.MassTable[1]= 1;
  All.MassTable[2]= 1;
  All.MassTable[3]= 1;
  All.MassTable[4]= 0;
  for(i=0, massflag=0;i<5;i++)
    {
      if(All.MassTable[i]==0 )
	massflag=1;
    }

  fprintf(stderr, "\nN_sph: %d\nN_halo: %d\nN_disk: %d\nN_bulge: %d\nN_stars: %d\n",
	  All.TotN_gas, All.TotN_halo, All.TotN_disk, All.TotN_bulge, 
	  All.TotN_stars);
       
    NumPart = All.TotNumPart =    All.TotN_gas  + All.TotN_halo 
      + All.TotN_disk + All.TotN_bulge + All.TotN_stars;

    All.MaxPart =  (int)All.PartAllocFactor *  All.TotNumPart;    
    /* sets the maximum number of particles that may 
       reside on a processor */
    All.MaxPartSph=  (int)All.PartAllocFactor * All.TotN_gas;

    fprintf(stderr, "Numpart=%d\n", NumPart);

    allocate_memory();
      
    fprintf(stderr, "allocated memory\n");

    if(massflag)
      {
	for(type=0, count=1; type<5; type++)
          {
            if(All.MassTable[type]==0 )
	      {
		for_all_daughters(dyn, b, bi) 
                  {
		    P[count++].Mass =  (float) bi->get_mass() ;
		  } } } }

    count=1;
    for_all_daughters(dyn, b, bi) 
      {
	du0 = bi->get_pos();

	for(k=0;k<3;k++) 
	  { P[count].PosPred[k]= (float) du0[k]; }
	count++;      
      }

    count=1;
    for_all_daughters(dyn, b, bi) 
      {
	du0 = bi->get_vel();

	for(k=0;k<3;k++) 
	  { P[count].VelPred[k]= (float) du0[k]; }
	count++;      
      }

    count=1;
    for_all_daughters(dyn, b, bi) 
      {
	P[count].ID=count++; 
      }
      
    fprintf(stderr,"done with reading.\n"); fflush(stderr);

      
    /* set the particle types */
    for(type=4, pc=1; type<5; type++)
      for(i=0; i<n; i++)
	P[pc++].Type = type;

    /* set temperature if desired */

    fprintf(stderr,"Baryonic particles        :  %d\n", N_gas);
    fprintf(stderr,"Collisionless particles   :  %d\n", NumPart-N_gas);
    fprintf(stderr,"                          ----------\n");
    fprintf(stderr,"Total number of particles :  %d\n\n", NumPart);
  
    savepositions(1);
  
}


/* This routine allocates memory for 
 * particle storage, both the collisionless and the SPH particles.
 * The memory for the ordered binary tree of the timeline
 * is also allocated.
 */
void allocate_memory(void)
{
  int bytes,bytes_tot=0;

  if(All.MaxPart>0)
    {
      if(!(P_data= (struct particle_data*) malloc(bytes=All.MaxPart*sizeof(struct particle_data))))
        {
          fprintf(stderr, "failed to allocate memory for `P_data' (%d bytes).\n",bytes);
          endrun(1);
        }
      bytes_tot+=bytes;

      P= P_data-1;   /* start with offset 1 */

      fprintf(stderr, "\nAllocated %g MByte for particle storage.\n\n",bytes_tot/(1024.0*1024.0));
    }
}



/* This function writes a snapshot of the particle ditribution to
 * one file using Gadget's default file format.
 * Each snapshot file contains a header first, then particle positions, 
 * velocities and ID's.
 * Then particle masses are written for those particle types with zero entry in
 * MassTable.
 * After that, first the internal energies u, and then the density is written
 * for the SPH particles.
 * Finally, if cooling is enabled, the mean molecular weight is written for the gas
 * particles. 
 */

void savepositions(int num)
{
  FILE *fd;
  char buf[100];
  float dummy[3];
  int i,k;
  int   blklen,masscount;
  double a3inv;
#ifdef COOLING
  double ne, nh0;
#endif

#define BLKLEN my_fwrite(&blklen, sizeof(blklen), 1, fd);


  if(All.ComovingIntegrationOn)
    a3inv=  1/(All.Time*All.Time*All.Time);
  else
    a3inv=  1.0;

  //  sprintf(buf,"%s","GDT_format_datafile");

  if((fd=stdout))
    {
      header1.npart[0]= header1.npartTotal[0]= All.TotN_gas;
      header1.npart[1]= header1.npartTotal[1]= All.TotN_halo;
      header1.npart[2]= header1.npartTotal[2]= All.TotN_disk;
      header1.npart[3]= header1.npartTotal[3]= All.TotN_bulge;
      header1.npart[4]= header1.npartTotal[4]= All.TotN_stars;
      header1.npart[5]= header1.npartTotal[5]= 0;

      for(i=0;i<6;i++)
        header1.mass[i]=0;

      for(i=0, masscount=0; i<5; i++)
        {
          header1.mass[i]= All.MassTable[i];
          if(All.MassTable[i]==0 && header1.npart[i]>0)
            masscount+= header1.npart[i];
        }

      header1.time= All.Time;

      if(All.ComovingIntegrationOn)
        header1.redshift=1.0/All.Time - 1.0;
      else
        header1.redshift=0;  

      
      header1.flag_sfr=0;
      header1.flag_feedback=0;
      header1.flag_cooling= 0;
#ifdef COOLING
      header1.flag_cooling= 1;
#endif
      header1.num_files= 1;
      header1.BoxSize= All.BoxSize;
      header1.Omega0=  All.Omega0;
      header1.OmegaLambda= All.OmegaLambda;
      header1.HubbleParam= All.HubbleParam;
      
      blklen=sizeof(header1);
      BLKLEN;
      my_fwrite(&header1, sizeof(header1), 1, fd);
      BLKLEN;


      blklen=NumPart*3*sizeof(float);

      BLKLEN;
      for(i=1;i<=NumPart;i++)
        {
          for(k=0;k<3;k++)
            dummy[k]=P[i].PosPred[k];
          my_fwrite(dummy,sizeof(float),3,fd);
        }
      BLKLEN;


      BLKLEN;
      for(i=1;i<=NumPart;i++)
        {
          for(k=0;k<3;k++)
            dummy[k]=P[i].VelPred[k];
          my_fwrite(dummy,sizeof(float),3,fd);
        }
      BLKLEN;
 

      blklen=NumPart*sizeof(int);
      BLKLEN;
      for(i=1;i<=NumPart;i++)
        {
          my_fwrite(&P[i].ID,sizeof(int),1,fd);
        }
      BLKLEN;

      blklen=masscount*sizeof(float);
      if(masscount)
        BLKLEN;
      for(i=1;i<=NumPart;i++)
        {
          dummy[0]= P[i].Mass;
          if(All.MassTable[P[i].Type]==0)
            my_fwrite(dummy,sizeof(float),1,fd);
        }
      if(masscount)
        BLKLEN;

      if(N_gas)
        {
          blklen=N_gas*sizeof(float);
          BLKLEN;
          for(i=1;i<=N_gas;i++)
            {
              dummy[0]=SphP[i].EgySpecPred;
              my_fwrite(dummy,sizeof(float),1,fd);
            }
          BLKLEN;


          blklen=N_gas*sizeof(float);  /* added density  */
          BLKLEN;
          for(i=1;i<=N_gas;i++)
            {
              dummy[0]=SphP[i].DensityPred;
              my_fwrite(dummy,sizeof(float),1,fd);
            }
          BLKLEN;

#ifdef COOLING
          blklen=N_gas*sizeof(float);  /* electron abundance */
          BLKLEN;
          for(i=1;i<=N_gas;i++)
            {
              dummy[0]= SphP[i].Ne;
              my_fwrite(dummy,sizeof(float),1,fd);
            }
          BLKLEN;


          blklen=N_gas*sizeof(float);  /* neutral hydrogen */
          BLKLEN;
          for(i=1;i<=N_gas;i++)
            {
              ne= SphP[i].Ne;

              AbundanceRatios(SphP[i].EgySpecPred, SphP[i].DensityPred*a3inv,
                              &ne, &nh0);
              dummy[0]= nh0;
              my_fwrite(dummy,sizeof(float),1,fd);
            }
          BLKLEN;
#endif
        }
      fclose(fd);
    }
  else
    {
      fprintf(stderr,"Error. Can't write in file '%s'\n", buf);
      endrun(10);
    }
}


/* This catches I/O errors occuring for my_my_fwrite(). In this case we better stop.
 */
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t nwritten;

  if((nwritten=fwrite(ptr, size, nmemb, stream))!=nmemb)
    {
      fprintf(stderr, "I/O error (fwrite) on has occured.\n");
      fflush(stderr);
      endrun(777);
    }
  return nwritten;
}


/* This catches I/O errors occuring for fread(). In this case we better stop.
 */
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t nread;

  if((nread=fread(ptr, size, nmemb, stream))!=nmemb)
    {
      fprintf(stderr, "I/O error (fread) has occured.\n");
      fflush(stderr);
      endrun(778);
    }
  return nread;
}


void endrun(int ierr)
{
  if(ierr)
    {
      fprintf(stderr,"endrun called with an error level of %d\n\n\n", ierr);
      exit(1);
    }
  exit(0);
}


void main(int argc, char ** argv)
{
  check_help();

  extern char *poptarg;
  int c;
  char* param_string = "iw";

  bool i_flag = false;
  bool w_flag = false;

  while ((c = pgetopt(argc, argv, param_string,
		      "$Revision$", _SRC_)) != -1)
    switch(c) {

    case 'i': i_flag = true;
      break;
    case 'w': w_flag = true;
      break;
    case '?': params_to_usage(cerr, argv[0], param_string);
      exit(1);
    }

  dyn *root;

  if(w_flag) {
    write_gadget();
  }
  else {
    root = read_gadget(i_flag);
    root->log_history(argc, argv);
    put_node(root, cout);
  }

}

#endif

/* end of: readgadget.C */

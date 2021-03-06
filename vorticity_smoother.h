


// namespace extension
namespace oomph
{

//===============================================
/// Overloaded element that allows projection of
/// vorticity.
//===============================================
template<class ELEMENT>
class VorticitySmootherElement : public virtual ELEMENT
{

public:


 /// Constructor
 VorticitySmootherElement()
  {
   // hierher only coded up for 2D at the moment. Check.

   // Index of smoothed vorticity 
   Smoothed_vorticity_index=3;

   // Number of values per field: u,v,p,omega,d/dx,d/dy, 
   // d^2/dx^2,d^2/dxdy, d^2/dy^2, 
   // d^3/dx^3, d^3/dx^2dy, d^3/dxdy^2,d^3,dy^3
   // du/dx, du/dy, dv/dx, dv/dy
   Number_of_values_per_field=17; // hierher 13; 

   // Pointer to fct that specifies exact vorticity and
   // derivs (for validation).
   Exact_vorticity_fct_pt=0;
  } 

 /// Typedef for pointer to function that specifies the exact
 /// vorticity and derivs (for validation)
 typedef void (*ExactVorticityFctPt)(const Vector<double>& x, 
                                     double& vort,
                                     Vector<double>& dvort_dx,
                                     Vector<double>& dvort_dxdy,
                                     Vector<double>& dvort_dxdxdy,
                                     Vector<double>& dveloc_dx);
  
 /// Access function: Pointer to fct that specifies exact vorticity and
 /// derivs (for validation).
 ExactVorticityFctPt& exact_vorticity_fct_pt() 
  {return Exact_vorticity_fct_pt;}

 /// Access function: Pointer to fct that specifies exact vorticity and
 /// derivs (for validation). const version
 ExactVorticityFctPt exact_vorticity_fct_pt() const 
 {return Exact_vorticity_fct_pt;}


 /// Index of smoothed vorticity -- followed by derivatives
 unsigned smoothed_vorticity_index() const
 {
  return Smoothed_vorticity_index;
 }
 

 /// \short Number of values required at local node n. In order to simplify
 /// matters, we allocate storage for pressure variables at all the nodes
 /// and then pin those that are not used. 
 unsigned required_nvalue(const unsigned &n) const 
 {return Number_of_values_per_field;} 


 /// \short Number of continuously interpolated values:
 unsigned ncont_interpolated_values() const 
 {return Number_of_values_per_field;} 

/// \short Get the function value u in Vector.
/// Note: Given the generality of the interface (this function
/// is usually called from black-box documentation or interpolation routines),
/// the values Vector sets its own size in here.
 void get_interpolated_values(const Vector<double>&s,  Vector<double>& values)
  {
   unsigned t=0;
   get_interpolated_values(t,s,values);
  }
 
/// \short Get the function value u in Vector.
/// Note: Given the generality of the interface (this function
/// is usually called from black-box documentation or interpolation routines),
/// the values Vector sets its own size in here.
 void get_interpolated_values(const unsigned& t, const Vector<double>&s, 
                              Vector<double>& values)
  {
   unsigned DIM=2;

   // Set size of Vector
   values.resize(Number_of_values_per_field); 
   
   // Initialise
   for(unsigned i=0;i<DIM+1;i++) {values[i]=0.0;} 
   
   //Find out how many nodes there are
   unsigned n_node = this->nnode();
   
   // Shape functions
   Shape psif(n_node);
   this->shape(s,psif);
   
   //Calculate velocities: values[0],...
   for(unsigned i=0;i<DIM;i++) 
    {
     //Get the index at which the i-th velocity is stored
     unsigned u_nodal_index = this->u_index_nst(i);
     for(unsigned l=0;l<n_node;l++) 
      {
       values[i] += this->nodal_value(t,l,u_nodal_index)*psif[l];
      } 
    }
   
   //Calculate pressure: values[DIM] 
   //(no history is carried in the pressure)
   values[DIM] = this->interpolated_p_nst(s);


   // No need to interpolate these onto the new mesh
   /* // Vorticity and its derivs */
   /* for (unsigned i=0;i<6;i++)  */
   /*  { */
   /*   for(unsigned l=0;l<n_node;l++)  */
   /*    { */
   /*     values[3+i]+=this->nodal_value(t,l,Smoothed_vorticity_index+i)*psif[l]; */
   /*    }  */
   /*  } */

  }

 /// Pin all smoothed vorticity quantities
 void pin_smoothed_vorticity()
  {
   unsigned nnod=this->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt=this->node_pt(j);
     for (unsigned i=Smoothed_vorticity_index;
          i<Number_of_values_per_field;i++)
      {
       nod_pt->pin(i);
      }
    }
  }


 /// \short Output exact veloc, vorticity, derivs and indicator
 /// based on functions specified by two function pointers
 void output_analytical_veloc_and_vorticity(std::ostream &outfile, 
                                            const unsigned &nplot)
 {
   //Vector of local coordinates
   Vector<double> s(2);
   
   // Shape functions
   unsigned n_node = this->nnode();   
   Shape psif(n_node);

   // Tecplot header info
   outfile << this->tecplot_zone_string(nplot);
   
   // Loop over plot points
   unsigned num_plot_points=this->nplot_points(nplot);
   for (unsigned iplot=0;iplot<num_plot_points;iplot++)
    {
     // Get local coordinates of plot point
     this->get_s_plot(iplot,nplot,s);

     // Coordinates
     Vector<double> x(2);
     for(unsigned i=0;i<2;i++) 
      {
       x[i]=this->interpolated_x(s,i);
       outfile << x[i] << " ";
      }

     // Fake veloc and  pressure
     outfile << "0.0 0.0 0.0 ";

     // Get vorticity and its derivatives
     double vort=0.0;
     Vector<double> dvort_dx(2);
     Vector<double> dvort_dxdy(3);
     Vector<double> dvort_dxdxdy(4);
     Vector<double> dveloc_dx(4);
     Exact_vorticity_fct_pt(x,
                            vort,
                            dvort_dx,
                            dvort_dxdy,
                            dvort_dxdxdy,
                            dveloc_dx);
     
     // Smoothed vorticity
     outfile << vort << " ";

     // Smoothed vorticity derivatives (d/dx, d/dy, d^2/dx^2, d^2/dxdy, d^2/dy^2
     // d^3/dx^3, d^3/dx^2dy, d^3/dxdy^2, d^3/dy^3,
     // du/dx, du/dy, dv/dx, dv/dy
     outfile << dvort_dx[0] << " "
             << dvort_dx[1] << " "
             << dvort_dxdy[0] << " "
             << dvort_dxdy[1] << " "
             << dvort_dxdy[2] << " "
             << dvort_dxdxdy[0] << " "
             << dvort_dxdxdy[1] << " "
             << dvort_dxdxdy[2] << " "
             << dvort_dxdxdy[3] << " "
             << dveloc_dx[0] << " "
             << dveloc_dx[1] << " "
             << dveloc_dx[2] << " "
             << dveloc_dx[3] << " ";
     
     outfile << std::endl;   
    }
   
   // Write tecplot footer (e.g. FE connectivity lists)
   this->write_tecplot_zone_footer(outfile,nplot);
  }




 /// \short Output veloc, smoothed vorticity and derivatives
 void output_smoothed_vorticity(std::ostream &outfile, 
                                const unsigned &nplot)
  {
   //Vector of local coordinates
   Vector<double> s(2);
   
   // Shape functions
   unsigned n_node = this->nnode();   
   Shape psif(n_node);

   // Tecplot header info
   outfile << this->tecplot_zone_string(nplot);
   
   // Loop over plot points
   unsigned num_plot_points=this->nplot_points(nplot);
   for (unsigned iplot=0;iplot<num_plot_points;iplot++)
    {
     // Get local coordinates of plot point
     this->get_s_plot(iplot,nplot,s);

     // Get vorticity and its derivatives (reconstructed)
     double vort=0.0;
     Vector<double> veloc(2);
     Vector<double> dvort_dx(2);
     Vector<double> dvort_dxdy(3);
     Vector<double> dvort_dxdxdy(4); 
     Vector<double> dveloc_dx(4);
     vorticity_and_its_derivs(s, 
                              veloc,
                              vort,
                              dvort_dx,
                              dvort_dxdy,
                              dvort_dxdxdy,
                              dveloc_dx);
     // Coordinates
     Vector<double> x(2);
     for(unsigned i=0;i<2;i++) 
      {
       x[i]=this->interpolated_x(s,i);
       outfile << x[i] << " ";
      }
     
     // Veloc
     outfile << veloc[0] << " "
             << veloc[1] << " ";
     
     // Smoothed vorticity
     outfile << vort << " ";

     // Smoothed vorticity derivatives (d/dx, d/dy, d^2/dx^2, d^2/dxdy, d^2/dy^2
     // d^3/dx^3, d^3/dx^2dy, d^3/dxdy^2, d^3/dy^3
     // du/dx, du/dy, dv/dx, dv/dy
     outfile << dvort_dx[0] << " "
             << dvort_dx[1] << " "
             << dvort_dxdy[0] << " "
             << dvort_dxdy[1] << " "
             << dvort_dxdy[2] << " "
             << dvort_dxdxdy[0] << " "
             << dvort_dxdxdy[1] << " "
             << dvort_dxdxdy[2] << " "
             << dvort_dxdxdy[3] << " "
             << dveloc_dx[0] << " "
             << dveloc_dx[1] << " "
             << dveloc_dx[2] << " "
             << dveloc_dx[3] << " ";

     outfile << std::endl;   
    }
   
   // Write tecplot footer (e.g. FE connectivity lists)
   this->write_tecplot_zone_footer(outfile,nplot);
  }







 /// \short Overloaded output fct: Output veloc, pressure, 
 /// smoothed vorticity 
 void output(std::ostream &outfile, const unsigned &nplot)
  {
   //Vector of local coordinates
   Vector<double> s(2);
   
   // Shape functions
   unsigned n_node = this->nnode();
   Shape psif(n_node);

   // Tecplot header info
   outfile << this->tecplot_zone_string(nplot);
   
   // Loop over plot points
   unsigned num_plot_points=this->nplot_points(nplot);
   for (unsigned iplot=0;iplot<num_plot_points;iplot++)
    {
     // Get local coordinates of plot point
     this->get_s_plot(iplot,nplot,s);
     
     // Get shape fct
     this->shape(s,psif);

     // Coordinates
     Vector<double> x(2);
     for(unsigned i=0;i<2;i++)
      {
       x[i]=this->interpolated_x(s,i);
       outfile << x[i] << " ";
      }
     
     // Veloc
     Vector<double> veloc(2);
     this->interpolated_u_nst(s,veloc);
     outfile << veloc[0] << " "
             << veloc[1] << " ";

     // Pressure
     outfile << this->interpolated_p_nst(s) << " ";

     // Smoothed vorticity
     double smoothed_vort=0.0;
     for(unsigned l=0;l<n_node;l++)
      {
       smoothed_vort += this->nodal_value(l,Smoothed_vorticity_index)*psif[l];
      }
     outfile << smoothed_vort << " ";

     // Smoothed vorticity derivatives (d/dx, d/dy, d^2/dx^2, d^2/dxdy, d^2/dy^2
     // d^3/dx^3, d^3/dx^2dy, d^3/dxdy^2, d^3/dy^3,
     // du/dx, du/dy, dv/dx, dv/dy
     for (unsigned i=1;i<14;i++) // hierher was 10
      {
       double smoothed_vort_deriv=0.0;
       for(unsigned l=0;l<n_node;l++)
        {
         smoothed_vort_deriv+=
          this->nodal_value(l,Smoothed_vorticity_index+i)*psif[l];
        }
       outfile << smoothed_vort_deriv << " ";
      }

     outfile << std::endl;
    }
   
   // Write tecplot footer (e.g. FE connectivity lists)
   this->write_tecplot_zone_footer(outfile,nplot);
  }


 /// Get raw derivative of velocity
 void get_raw_velocity_deriv(const Vector<double>& s,
                             Vector<double>& dveloc_dx) const
 {
  //Find out how many nodes there are
  unsigned n_node = this->nnode();
  
  //Set up memory for the shape functions
  Shape psif(n_node);
  DShape dpsifdx(n_node,2);
  
  //Call the derivatives of the shape and test functions
  this->dshape_eulerian(s,psif,dpsifdx);
  
  //Initialise to zero
  for(unsigned j=0;j<4;j++)
   {
    dveloc_dx[j] = 0.0;
   }
  
  // Loop over nodes
  for(unsigned l=0;l<n_node;l++) 
   {
    //Loop over derivative directions
    for(unsigned j=0;j<2;j++)
     {                               
      dveloc_dx[j]   += this->nodal_value(l,0)*dpsifdx(l,j);            
      dveloc_dx[j+2] += this->nodal_value(l,1)*dpsifdx(l,j);
      }
    }
  }

 /// Get raw derivative of smoothed vorticity
 void get_raw_vorticity_deriv(const Vector<double>& s,
                          Vector<double>& dvorticity_dx) const
  {
   //Find out how many nodes there are
   unsigned n_node = this->nnode();
   
   //Set up memory for the shape functions
   Shape psif(n_node);
   DShape dpsifdx(n_node,2);

   //Call the derivatives of the shape and test functions
   this->dshape_eulerian(s,psif,dpsifdx);

   //Initialise to zero
   for(unsigned j=0;j<2;j++)
    {
     dvorticity_dx[j] = 0.0;
    }
   
   // Loop over nodes
   for(unsigned l=0;l<n_node;l++) 
    {
     //Loop over derivative directions
     for(unsigned j=0;j<2;j++)
      {                               
       dvorticity_dx[j] += this->nodal_value(l,Smoothed_vorticity_index)*
        dpsifdx(l,j);
      }
    }
  }

 /// Get raw derivative of smoothed derivative vorticity
 void get_raw_vorticity_second_deriv(const Vector<double>& s,
                                     Vector<double>& dvorticity_dxdy) const
 {
  //Find out how many nodes there are
  unsigned n_node = this->nnode();
  
  //Set up memory for the shape functions
  Shape psif(n_node);
  DShape dpsifdx(n_node,2);
  
  //Call the derivatives of the shape and test functions
  this->dshape_eulerian(s,psif,dpsifdx);
  
  //Initialise to zero
  for(unsigned j=0;j<3;j++)
   {
    dvorticity_dxdy[j] = 0.0;
   }
  
  // Loop over nodes
  for(unsigned l=0;l<n_node;l++) 
   {
    //Loop over derivative directions to obtain xx and xy derivatives
    for(unsigned j=0;j<2;j++)
     {                               
      dvorticity_dxdy[j] += this->nodal_value(l,Smoothed_vorticity_index+1)*
       dpsifdx(l,j);
     }
    //Calcutation of yy derivative
    dvorticity_dxdy[2] += this->nodal_value(l,Smoothed_vorticity_index+2)*
     dpsifdx(l,1);
   }
 }
 

 /// Get raw derivative of smoothed derivative vorticity
 /// [0]: d^3/dx^3, [0]: d^3/dx^2dy, [0]: d^3/dxdy^2, [0]: d^3/dy^3, 
 void get_raw_vorticity_third_deriv(const Vector<double>& s,
                                    Vector<double>& dvorticity_dxdy) const
 {

  //Find out how many nodes there are
  unsigned n_node = this->nnode();
  
  //Set up memory for the shape functions
  Shape psif(n_node);
  DShape dpsifdx(n_node,2);
  
  //Call the derivatives of the shape and test functions
  this->dshape_eulerian(s,psif,dpsifdx);
  
  //Initialise to zero
  for(unsigned j=0;j<4;j++)
   {
    dvorticity_dxdy[j] = 0.0;
   }
  
  // Loop over nodes
  for(unsigned l=0;l<n_node;l++) 
   {
    // d^3/dx^3 = d/dx \overline{d^2/dx^2} 
    dvorticity_dxdy[0] += this->nodal_value(l,Smoothed_vorticity_index+3)*
     dpsifdx(l,0);
     
    // d^3/dx^2dy = d/dx \overline{d^2/dxdy} 
    dvorticity_dxdy[1] += this->nodal_value(l,Smoothed_vorticity_index+4)*
     dpsifdx(l,0);

    // d^3/dxdy^2 = d/dy \overline{d^2/dxdy} 
    dvorticity_dxdy[2] += this->nodal_value(l,Smoothed_vorticity_index+4)*
     dpsifdx(l,1);

    // d^3/dy^3 = d/dy \overline{d^2/dy^2} 
    dvorticity_dxdy[3] += this->nodal_value(l,Smoothed_vorticity_index+5)*
     dpsifdx(l,1);

   }
 }
 

 /// \short Compute the element's contribution to the (squared) L2 norm
 /// of the difference between exact and smoothed vorticity. i=0: do 
 /// vorticity itself; i>0: derivs 
 double vorticity_error_squared(const unsigned& i)
 {  
  double norm_squared=0.0;
  
  //Find out how many nodes there are
  unsigned n_node = this->nnode();
  
  //Set up memory for the shape functions
  Shape psif(n_node);
  DShape dpsifdx(n_node,2);
  
  //Number of integration points
  unsigned n_intpt = this->integral_pt()->nweight();
  
  //Set the Vector to hold local coordinates
  Vector<double> s(2);
  Vector<double> x(2);
  
  //Loop over the integration points
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
   {
    //Assign values of s
    for(unsigned ii=0;ii<2;ii++)
     {
      x[ii]=0.0;
      s[ii] = this->integral_pt()->knot(ipt,ii);
     }
    
    //Get the integral weight
    double w = this->integral_pt()->weight(ipt);
    
    //Call the derivatives of the shape and test functions
    double J = this->dshape_eulerian(s,psif,dpsifdx);
    
    //Premultiply the weights and the Jacobian
    double W = w*J;
    
    // Smoothed vorticity
    double smoothed_vort=0.0;
    for(unsigned l=0;l<n_node;l++) 
     {
      Node* nod_pt=this->node_pt(l);
      smoothed_vort += this->nodal_value(l,Smoothed_vorticity_index+i)*psif[l];
      for (unsigned ii=0;ii<2;ii++)
       {
        x[ii]+=nod_pt->x(ii)*psif[l];
       }
     }   
    
     // Synthetic quantities
     double synth_vort=0.0;
     Vector<double> synth_dvort_dx(2,0.0);
     Vector<double> synth_dvort_dxdy(3,0.0);
     Vector<double> synth_dvort_dxdxdy(4,0.0); 
     Vector<double> synth_dveloc_dx(4,0.0);

     if (Exact_vorticity_fct_pt!=0) 
      {
       Exact_vorticity_fct_pt(x,synth_vort,synth_dvort_dx,
                              synth_dvort_dxdy,synth_dvort_dxdxdy,
                              synth_dveloc_dx);
      }
    double synth_quantity=synth_vort;
    if ((i==1)||(i==2))
     {
      synth_quantity=synth_dvort_dx[i-1];
     }    
    if ((i==3)||(i==4)||(i==5))
     {
      synth_quantity=synth_dvort_dxdy[i-3];
     }
    if ((i==6)||(i==7)||(i==8)||(i==9))
     {
      synth_quantity=synth_dvort_dxdxdy[i-6];
     }
    if ((i==10)||(i==11)||(i==12)||(i==13))
     {
      synth_quantity=synth_dveloc_dx[i-10];
     }
    // Add squared difference
    norm_squared+=pow(smoothed_vort-synth_quantity,2)*W;
   }
  
  return norm_squared;
 }


 /// \short Compute smoothed vorticity and its derivatives
 void vorticity_and_its_derivs(const Vector<double>& s, 
                               Vector<double>& veloc, 
                               double& vort,
                               Vector<double>& dvort_dx,
                               Vector<double>& dvort_dxdy,
                               Vector<double>& dvort_dxdxdy,
                               Vector<double>& dveloc_dx) 
 {
  // Shape functions
  unsigned n_node = this->nnode();   
  Shape psif(n_node);
  this->shape(s,psif);
  
  // Smoothed vorticity
  vort=0.0;
  veloc[0]=0.0;
  veloc[1]=0.0;
  for(unsigned l=0;l<n_node;l++) 
   {
    vort += this->nodal_value(l,Smoothed_vorticity_index)*psif[l];
    veloc[0]+=this->nodal_value(l,0)*psif[l];
    veloc[1]+=this->nodal_value(l,1)*psif[l];
   }   
  
  // Smoothed vorticity derivatives (d/dx, d/dy, d^2/dx^2, d^2/dxdy, d^2/dy^2, 
  // d^3/dx^3, d^3/dx^2dy, d^3/dxdy^2,d^3,dy^3
  // du/dx, du/dy, dv/dx, dv/dy
  for (unsigned i=1;i<14;i++) // hierher used to be 10
   {
    double smoothed_vort_deriv=0.0;
    for(unsigned l=0;l<n_node;l++) 
     {
      smoothed_vort_deriv+=
       this->nodal_value(l,Smoothed_vorticity_index+i)*psif[l];
     } 
    switch (i)
     {
     case 1:
      dvort_dx[0]=smoothed_vort_deriv;
      break;
     case 2:
      dvort_dx[1]=smoothed_vort_deriv;
      break;
     case 3:
      dvort_dxdy[0]=smoothed_vort_deriv;
      break;
     case 4:
      dvort_dxdy[1]=smoothed_vort_deriv;
      break;
     case 5:
      dvort_dxdy[2]=smoothed_vort_deriv;
      break;
     case 6:
      dvort_dxdxdy[0]=smoothed_vort_deriv;
      break;
     case 7:
      dvort_dxdxdy[1]=smoothed_vort_deriv;
      break;
     case 8:
      dvort_dxdxdy[2]=smoothed_vort_deriv;
      break;
     case 9:
      dvort_dxdxdy[3]=smoothed_vort_deriv;
      break;
     case 10:
      dveloc_dx[0]=smoothed_vort_deriv;
      break;
     case 11:
      dveloc_dx[1]=smoothed_vort_deriv;
      break;
     case 12:
      dveloc_dx[2]=smoothed_vort_deriv;
      break;
     case 13:
      dveloc_dx[3]=smoothed_vort_deriv;
      break;
     default:
      oomph_info << "never get here\n";
      abort();
      break;
     }
   }
 }

  private:

 /// Index of smoothed vorticity -- followed by derivatives
 unsigned Smoothed_vorticity_index;

 /// \short Number of values per field: u,v,p,omega,d/dx,d/dy, 
 /// d^2/dx^2,d^2/dxdy, d^2/dy^2,
 // du/dx, du/dy, dv/dx, dv/dy
 unsigned Number_of_values_per_field;

 /// Pointer to fct that specifies exact vorticity and
 /// derivs (for validation).
 ExactVorticityFctPt Exact_vorticity_fct_pt;


};


} // end namespace extension






//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////




//========================================================
/// Smoother for vorticity in 2D
//========================================================
template<class ELEMENT>
class VorticitySmoother 
{
   public:
 
 
 /// Constructor: Set order of recovery shape functions
 VorticitySmoother(const unsigned& recovery_order) : 
  Recovery_order(recovery_order)
  {}
 
  /// Broken copy constructor
 VorticitySmoother(const VorticitySmoother&) 
  { 
   BrokenCopy::broken_copy("VorticitySmoother");
  } 
 
 /// Broken assignment operator
 void operator=(const VorticitySmoother&) 
  {
   BrokenCopy::broken_assign("VorticitySmoother");
  }
 
 /// Empty virtual destructor
 virtual ~VorticitySmoother(){}
 
 /// Access function for order of recovery polynomials
 unsigned& recovery_order() {return Recovery_order;}
 
 /// Recovery shape functions as functions of the global, Eulerian
 /// coordinate x of dimension dim.
 /// The recovery shape functions are  complete polynomials of 
 /// the order specified by Recovery_order.
 void shape_rec(const Vector<double>& x,
                Vector<double>& psi_r)
  {
   std::ostringstream error_stream;
   
   /// Find order of recovery shape functions
   switch(recovery_order())
    {
    case 1:
     
     // Complete linear polynomial in 2D:
     psi_r[0]=1.0;
     psi_r[1]=x[0];
     psi_r[2]=x[1];
     break;
     
    case 2:
     
     // Complete quadratic polynomial in 2D:
     psi_r[0]=1.0;
     psi_r[1]=x[0];
     psi_r[2]=x[1];
     psi_r[3]=x[0]*x[0];
     psi_r[4]=x[0]*x[1];
     psi_r[5]=x[1]*x[1];
     break;
     
    case 3:
     
     // Complete cubic polynomial in 2D:
     psi_r[0]=1.0;
     psi_r[1]=x[0];
     psi_r[2]=x[1];
     psi_r[3]=x[0]*x[0];
     psi_r[4]=x[0]*x[1];
     psi_r[5]=x[1]*x[1];
     psi_r[6]=x[0]*x[0]*x[0];
     psi_r[7]=x[0]*x[0]*x[1];
     psi_r[8]=x[0]*x[1]*x[1];
     psi_r[9]=x[1]*x[1]*x[1];
     break;
     
    default:
     
     error_stream 
      << "Recovery shape functions for recovery order " 
      << recovery_order() << " haven't yet been implemented for 2D" 
      << std::endl;
     
     throw OomphLibError(error_stream.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
  }



/// Integation scheme associated with the recovery shape functions
/// must be of sufficiently high order to integrate the mass matrix
/// associated with the recovery shape functions. The  argument
/// is the dimension of the elements.
/// The integration is performed locally over the elements, so the 
/// integration scheme does depend on the geometry of the element.
/// The type of element is specified by the boolean which is
/// true if elements in the patch are QElements and false if they are 
/// TElements (will need change if we ever have other element types)
 Integral* integral_rec(const bool &is_q_mesh)
 {
  std::ostringstream error_stream; 
  
  // 2D:
  
  /// Find order of recovery shape functions
  switch(recovery_order())
   {
   case 1:
    
    // Complete linear polynomial in 2D:
    if(is_q_mesh) {return(new Gauss<2,2>);}
    else {return(new TGauss<2,2>);}
    break;
    
   case 2:
   
    // Complete quadratic polynomial in 2D:
    if(is_q_mesh) {return(new Gauss<2,3>);}
    else {return(new TGauss<2,3>);}
    break;
   
   case 3:
   
    // Complete cubic polynomial in 2D:
    if(is_q_mesh) {return(new Gauss<2,4>);}
    else {return(new TGauss<2,4>);}
    break;
   
   default:
   
    error_stream 
     << "Recovery shape functions for recovery order " 
     << recovery_order() << " haven't yet been implemented for 2D" 
     << std::endl;
   
    throw OomphLibError(error_stream.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
 
 
  //Dummy return (never get here)
  return 0;
 }
 
 
 /// \short Setup patches: For each vertex node pointed to by nod_pt,
 /// adjacent_elements_pt[nod_pt] contains the pointer to the vector that 
 /// contains the pointers to the elements that the node is part of.
 /// Also returns a Vector of vertex nodes for use in get_element_errors.
 void setup_patches(Mesh*& mesh_pt,
                    std::map<Node*,Vector<ELEMENT*>*>& adjacent_elements_pt,
                    Vector<Node*>& vertex_node_pt)
 {
   
  // Clear: hierher should we do this in Z2 as well?
  adjacent_elements_pt.clear();

  // Auxiliary map that contains element-adjacency for ALL nodes
  std::map<Node*,Vector<ELEMENT*>*> aux_adjacent_elements_pt;
   
#ifdef PARANOID
  // Check if all elements request the same recovery order
  unsigned ndisagree=0;
#endif
   
  // Loop over all elements to setup adjacency for all nodes.
  // Need to do this because midside nodes can be corner nodes for 
  // adjacent smaller elements! Admittedly, the inclusion of interior 
  // nodes is wasteful...
  unsigned nelem=mesh_pt->nelement();
  for (unsigned e=0;e<nelem;e++)
   {
    ELEMENT* el_pt=
     dynamic_cast<ELEMENT*>(mesh_pt->element_pt(e));
     
#ifdef PARANOID
    // Check if all elements request the same recovery order
    if (el_pt->nrecovery_order()!=Recovery_order){ndisagree++;}
#endif
     
    // Loop all nodes in element
    unsigned nnod=el_pt->nnode();
    for (unsigned n=0;n<nnod;n++)
     {
      Node* nod_pt=el_pt->node_pt(n);

      // Has this node been considered before?
      if (aux_adjacent_elements_pt[nod_pt]==0)
       {
        // Create Vector of pointers to its adjacent elements
        aux_adjacent_elements_pt[nod_pt]=new 
         Vector<ELEMENT*>;
       }
       
      // Add pointer to adjacent element
      (*aux_adjacent_elements_pt[nod_pt]).push_back(el_pt);       
     }
   } // end element loop
   
#ifdef PARANOID
  // Check if all elements request the same recovery order
  if (ndisagree!=0)
   {
    oomph_info 
     << "\n\n========================================================\n";
    oomph_info << "WARNING: " << std::endl;
    oomph_info << ndisagree << " out of " << mesh_pt->nelement() 
               << " elements\n";
    oomph_info << "have different preferences for the order of the recovery\n";
    oomph_info << "shape functions. We are using: Recovery_order=" 
               << Recovery_order << std::endl;
    oomph_info 
     << "========================================================\n\n";
   }
#endif
   
  //Loop over all elements, extract adjacency for corner nodes only
  nelem=mesh_pt->nelement();
  for (unsigned e=0;e<nelem;e++)
   {
    ELEMENT* el_pt=
     dynamic_cast<ELEMENT*>(mesh_pt->element_pt(e));
     
    // Loop over corner nodes
    unsigned n_node=el_pt->nvertex_node();
    for (unsigned n=0;n<n_node;n++)
     {
      Node* nod_pt=el_pt->vertex_node_pt(n);
       
      // Has this node been considered before?
      if (adjacent_elements_pt[nod_pt]==0)
       {

        // Add the node pointer to the vertex node container
        vertex_node_pt.push_back(nod_pt);
         
        // Create Vector of pointers to its adjacent elements
        adjacent_elements_pt[nod_pt]=new 
         Vector<ELEMENT*>;
         
        // Copy across:
        unsigned nel=(*aux_adjacent_elements_pt[nod_pt]).size();
        for (unsigned e=0;e<nel;e++)
         {
          (*adjacent_elements_pt[nod_pt]).push_back(
           (*aux_adjacent_elements_pt[nod_pt])[e]);
         }
       }
     }

   } // end of loop over elements
   

   // Cleanup
   for (typename std::map<Node*,Vector<ELEMENT*>*>::iterator it=
         aux_adjacent_elements_pt.begin();
        it!=aux_adjacent_elements_pt.end();it++)
    {
     delete it->second;
    }
   
 }
  


 /// Given the vector of elements that make up a patch, compute
 /// the vector of recovered vorticity coefficients and return
 /// a pointer to it. n_deriv indicates which derivative of the
 /// vorticity is supposed to be smoothed: 0: zeroth (i.e. vorticity
 /// itself; 1: d/dx; 2: d/dy; 3: d^2/dx^2; 4: d^2/dxdy 5: d^2/dy^2 
 /// 6: d^3/dx^3, 7: d^3/dx^2dy, 8: d^3/dxdy^2, 9: d^3/dy^3,
 /// 10: du/dx, 11: du/dy, 12: dv/dx, 13: dv/dy
 void get_recovered_vorticity_in_patch(
  const Vector<ELEMENT*>& patch_el_pt,
  const unsigned& num_recovery_terms, 
  Vector<double>*& recovered_vorticity_coefficient_pt,
  unsigned& n_deriv)
 {
  
  // Create/initialise matrix for linear system
  DenseDoubleMatrix recovery_mat(num_recovery_terms,num_recovery_terms,0.0);
  
  // Ceate/initialise vector for RHS
  Vector<double> rhs(num_recovery_terms,0.0);
  
  //Create a new integration scheme based on the recovery order
  //in the elements
  //Need to find the type of the element, default is to assume a quad
  bool is_q_mesh=true;
  
  //If we can dynamic cast to the TElementBase, then it's a triangle/tet
  //Note that I'm assuming that all elements are of the same geometry, but
  //if they weren't we could adapt...
  if(dynamic_cast<TElementBase*>(patch_el_pt[0])) {is_q_mesh=false;}
  
  Integral* const integ_pt = this->integral_rec(is_q_mesh);
  
  //Loop over all elements in patch to assemble linear system
  unsigned nelem=patch_el_pt.size();
  for (unsigned e=0;e<nelem;e++)
   {
    // Get pointer to element
    ELEMENT* const el_pt=patch_el_pt[e];
    
    // Create storage for the recovery shape function values 
    Vector<double> psi_r(num_recovery_terms);
    
    //Create vector to hold local coordinates
    Vector<double> s(2); 
    
    //Loop over the integration points
    unsigned n_intpt = integ_pt->nweight();   
    for(unsigned ipt=0;ipt<n_intpt;ipt++)
     {
      //Assign values of s, the local coordinate
      for(unsigned i=0;i<2;i++)
       {
        s[i] = integ_pt->knot(ipt,i);
       }
      
      //Get the integral weight
      double w = integ_pt->weight(ipt);
      
      //Jacobian of mapping
      double J = el_pt->J_eulerian(s);
     
     // Interpolate the global (Eulerian) coordinate
     Vector<double> x(2);
     el_pt->interpolated_x(s,x);
     
     // Premultiply the weights and the Jacobian
     // and the geometric jacobian weight (used in axisymmetric
     // and spherical coordinate systems) -- hierher really fct of x?
     // probably yes, actually).
     double W = w*J*(el_pt->geometric_jacobian(x));
     
     // Recovery shape functions at global (Eulerian) coordinate
     shape_rec(x,psi_r);
     
     // Get FE estimates for vorticity: 
     Vector<double> vorticity(1); 
     if(n_deriv==0)
      {
       el_pt->get_vorticity(s,vorticity);
      }

     // Get FE estimates for deriv of vorticity: 
     Vector<double> deriv_vorticity(2); 
     Vector<double> second_deriv_vorticity(3); 
     Vector<double> third_deriv_vorticity(4); 
     Vector<double> deriv_velocity(4); 
     if((n_deriv==1)||(n_deriv==2))
      {
       el_pt->get_raw_vorticity_deriv(s,deriv_vorticity);
      }
     // Get FE estimates for second deriv of vorticity: 
     else if((n_deriv==3)||(n_deriv==4)||(n_deriv==5))
      {
       el_pt->get_raw_vorticity_second_deriv(s,second_deriv_vorticity);
      }
     // Get FE estimates for third deriv of vorticity: 
     else if((n_deriv==6)||(n_deriv==7)||(n_deriv==8)||(n_deriv==9))
      {
       el_pt->get_raw_vorticity_third_deriv(s,third_deriv_vorticity);
      }
     // Get FE estimates for derivs of velocity
     else if((n_deriv==10)||(n_deriv==11)||(n_deriv==12)||(n_deriv==13))
      {
       el_pt->get_raw_velocity_deriv(s,deriv_velocity);
      }

     // Add elemental RHSs and recovery matrix to global versions
     //----------------------------------------------------------
     if(n_deriv==0)
      {
       // RHS 
       // Loop over the nodes for the test functions 
       for(unsigned l=0;l<num_recovery_terms;l++)
        {
         rhs[l] += vorticity[0]*psi_r[l]*W;
        }
      }
     // Get x-derivative of vorticity
     else if(n_deriv==1)
      {
       // RHS 
       // Loop over the nodes for the test functions 
       for(unsigned l=0;l<num_recovery_terms;l++)
        {
         rhs[l] += deriv_vorticity[0]*psi_r[l]*W;
        }
      }
     // Get y-derivative of vorticity
     else if(n_deriv==2)
      {
       // RHS 
       // Loop over the nodes for the test functions 
       for(unsigned l=0;l<num_recovery_terms;l++)
        {
         rhs[l] += deriv_vorticity[1]*psi_r[l]*W;
        }
      }
     // Get xx-derivative of vorticity
     else if(n_deriv==3)
      {
       // RHS 
       // Loop over the nodes for the test functions 
       for(unsigned l=0;l<num_recovery_terms;l++)
        {
         rhs[l] += second_deriv_vorticity[0]*psi_r[l]*W;
        }
      }
     // Get xy-derivative of vorticity
     else if(n_deriv==4)
      {
       // RHS 
       // Loop over the nodes for the test functions 
       for(unsigned l=0;l<num_recovery_terms;l++)
        {
         rhs[l] += second_deriv_vorticity[1]*psi_r[l]*W;
        }
      }
     // Get yy-derivative of vorticity
     else if(n_deriv==5)
      {
       // RHS 
       // Loop over the nodes for the test functions 
       for(unsigned l=0;l<num_recovery_terms;l++)
        {
         rhs[l] += second_deriv_vorticity[2]*psi_r[l]*W;
        }
      }
     // Get xxx-derivative of vorticity
     else if(n_deriv==6)
      {
       // RHS 
       // Loop over the nodes for the test functions 
       for(unsigned l=0;l<num_recovery_terms;l++)
        {
         rhs[l] += third_deriv_vorticity[0]*psi_r[l]*W;
        }
      }
     // Get xxy-derivative of vorticity
     else if(n_deriv==7)
      {
       // RHS 
       // Loop over the nodes for the test functions 
       for(unsigned l=0;l<num_recovery_terms;l++)
        {
         rhs[l] += third_deriv_vorticity[1]*psi_r[l]*W;
        }
      }
     // Get xxy-derivative of vorticity
     else if(n_deriv==8)
      {
       // RHS 
       // Loop over the nodes for the test functions 
       for(unsigned l=0;l<num_recovery_terms;l++)
        {
         rhs[l] += third_deriv_vorticity[2]*psi_r[l]*W;
        }
      }
     // Get xxy-derivative of vorticity
     else if(n_deriv==9)
      {
       // RHS 
       // Loop over the nodes for the test functions 
       for(unsigned l=0;l<num_recovery_terms;l++)
        {
         rhs[l] += third_deriv_vorticity[3]*psi_r[l]*W;
        }
      }
     // Get x-derivative of u-velocity
     else if(n_deriv==10)
      {
       // RHS 
       // Loop over the nodes for the test functions 
       for(unsigned l=0;l<num_recovery_terms;l++)
        {
         rhs[l] += deriv_velocity[0]*psi_r[l]*W;
        }
      }
     // Get y-derivative of u-velocity
     else if(n_deriv==11)
      {
       // RHS 
       // Loop over the nodes for the test functions 
       for(unsigned l=0;l<num_recovery_terms;l++)
        {
         rhs[l] += deriv_velocity[1]*psi_r[l]*W;
        }
      }
     // Get x-derivative of v-velocity
     else if(n_deriv==12)
      {
       // RHS 
       // Loop over the nodes for the test functions 
       for(unsigned l=0;l<num_recovery_terms;l++)
        {
         rhs[l] += deriv_velocity[2]*psi_r[l]*W;
        }
      }
     // Get y-derivative of v-velocity
     else if(n_deriv==13)
      {
       // RHS 
       // Loop over the nodes for the test functions 
       for(unsigned l=0;l<num_recovery_terms;l++)
        {
         rhs[l] += deriv_velocity[3]*psi_r[l]*W;
        }
      }
     else
      {
       oomph_info << "Never get here\n";
       abort();
      }


     // Loop over the nodes for the test functions 
     for(unsigned l=0;l<num_recovery_terms;l++)
      {
       //Loop over the nodes for the variables
       for(unsigned l2=0;l2<num_recovery_terms;l2++)
        { 
         //Add contribution to recovery matrix
         recovery_mat(l,l2)+=psi_r[l]*psi_r[l2]*W;
        }
      }      
    }   
  } // End of loop over elements that make up patch. 
 
 //Delete the integration scheme
 delete integ_pt;
 
 // Linear system is now assembled: Solve recovery system
 
 // LU decompose the recovery matrix
 recovery_mat.ludecompose();
 
 // Back-substitute
 recovery_mat.lubksub(rhs);
 
 // Now create a matrix to store the vorticity recovery coefficients.
 // Pointer to this matrix will be returned.
 recovered_vorticity_coefficient_pt =
  new Vector<double>(num_recovery_terms);
 
 // Copy coefficients
 for (unsigned icoeff=0;icoeff<num_recovery_terms;icoeff++)
  {
   (*recovered_vorticity_coefficient_pt)[icoeff]=rhs[icoeff]; 
  }
 
}

 // Get the recovery order
 unsigned nrecovery_order() const
  {
   switch(Recovery_order)
    {
    case 1:
     
     // Linear recovery shape functions
     //--------------------------------
     return 3; // 1, x, y
     break;
     
     
    case 2:
     
     // Quadratic recovery shape functions
     //-----------------------------------
     return 6; // 1, x, y, x^2, xy, y^2
     break;
     
    case 3:
     
     // Cubic recovery shape functions
     //--------------------------------
     return 10; // 1, x, y, x^2, xy, y^2, x^3, y^3, x^2 y, x y^2
     break;
     
    default:
     
     // Any other recovery order?
     //--------------------------
     std::ostringstream error_stream;
     error_stream 
      << "Wrong Recovery_order " << Recovery_order << std::endl;
     
     throw OomphLibError(error_stream.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
  }
 
 
 /// Recover vorticity from patches
 void recover_vorticity(Mesh* mesh_pt)
  {
   DocInfo doc_info;
   doc_info.disable_doc();
   recover_vorticity(mesh_pt,doc_info);
  }



 /// Recover vorticity from patches -- output intermediate steps
 /// to directory specified by DocInfo object
 void recover_vorticity(Mesh* mesh_pt, DocInfo& doc_info)
  {
   
   double t_start=TimingHelpers::timer();

   // Make patches
   //-------------
   std::map<Node*,Vector<ELEMENT*>*> adjacent_elements_pt;
   Vector<Node*> vertex_node_pt;
   setup_patches(mesh_pt,
                 adjacent_elements_pt,
                 vertex_node_pt);
   
   // Determine number of coefficients for expansion of recovered vorticity
   // Use complete polynomial of given order for recovery
   unsigned num_recovery_terms=nrecovery_order();
 
   // hierher get from element
   unsigned Smoothed_vorticity_index=3;  
   
   // Counter for averaging of recovered vorticity and its derivatives
   map<Node* , unsigned> count;

   // Loop over derivatives
   for (unsigned deriv=0;deriv<14;deriv++) // hierher used to be 10
    {

     // Storage for accumulated nodal vorticity (used to compute
     // nodal averages)
     map<Node*, double> averaged_recovered_vort;
     
     //Calculation of vorticity
     //------------------------
     
     // Do patch recovery
     // unsigned  counter=0;
     for (typename std::map<Node*,Vector<ELEMENT*>*>::iterator it=
           adjacent_elements_pt.begin();
          it!=adjacent_elements_pt.end();it++)
      {
       
       // Setup smoothed vorticity field for patches
       Vector<double>* recovered_vorticity_coefficient_pt;
       get_recovered_vorticity_in_patch(*(it->second),
                                        num_recovery_terms, 
                                        recovered_vorticity_coefficient_pt,
                                        deriv);
       
       // Now get the nodal average of the recovered vorticity
       // (nodes are generally part of multiple patches)
       
       //Loop over all elements to get recovered vorticity
       unsigned nelem=(*(it->second)).size();
       for (unsigned e=0;e<nelem;e++)
        {
         // Get pointer to element
         ELEMENT* const el_pt=(*(it->second))[e];
         
         // Get the number of nodes by element
         unsigned nnode_el=el_pt->nnode();
         for(unsigned j=0;j<nnode_el;j++)
          {
           //Get local coordinates
           Vector<double> s(2);
           Node* nod_pt=el_pt->node_pt(j);
           el_pt->local_coordinate_of_node(j,s);
           
           // Interpolate the global (Eulerian) coordinate
           Vector<double> x(2);
           el_pt->interpolated_x(s,x);
           
           // Recovery shape functions at global (Eulerian) coordinate
           Vector<double> psi_r(num_recovery_terms);
           shape_rec(x,psi_r);
           
           // Assemble recovered vorticity
           double recovered_vort=0.0;
           for (unsigned i=0;i<num_recovery_terms;i++)
            {
             recovered_vort+=(*recovered_vorticity_coefficient_pt)[i]*psi_r[i];
            }
           
           // Keep adding
           averaged_recovered_vort[nod_pt]+=recovered_vort;
           count[nod_pt]++;
          }
        }

       // Cleanup
       delete recovered_vorticity_coefficient_pt; 
       recovered_vorticity_coefficient_pt=0;
      }
     
     //Loop over all nodes to actually work out the average
     unsigned nnod=mesh_pt->nnode();
     for(unsigned j=0;j<nnod;j++)
      {
       Node* nod_pt=mesh_pt->node_pt(j);

       //Calculate the values of the smoothed vorticity 
       averaged_recovered_vort[nod_pt]/=count[nod_pt];
       
       //Assign smoothed vorticity to nodal values
       nod_pt->set_value(Smoothed_vorticity_index+deriv,
                         averaged_recovered_vort[nod_pt]);
      }
     
     
     // Start again
     count.clear();
     
     
    } // end of loop over derivatives

   // Cleanup
   for (typename std::map<Node*,Vector<ELEMENT*>*>::iterator it=
         adjacent_elements_pt.begin();
        it!=adjacent_elements_pt.end();it++)
    {
     delete it->second;
    }

   oomph_info << "Time for vorticity recovery: " 
              << TimingHelpers::timer()-t_start 
              << " sec " << std::endl;
  }

  private:

 /// Order of recovery polynomials
 unsigned Recovery_order;

};


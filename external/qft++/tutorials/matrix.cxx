/* Copyright 2008 Mike Williams (mwill@jlab.org)
 *
 * This file is part of qft++.
 *
 * qft++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * qft++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with qft++.  If not, see <http://www.gnu.org/licenses/>.
 */
//_____________________________________________________________________________
/** @file matrix.cxx
 *  @author Mike Williams
 *
 *  @brief Provides example usage of the matrix module
 *
 *  This tutorial gives a number of examples of how to perform common matrix
 *  operations using the qft++ matrix package.
 */
//_____________________________________________________________________________

#include <getopt.h>
#include "matrix.h"

using namespace std;
//_____________________________________________________________________________

/// Prints program usage to the screen
void PrintUsage();

/// Prints a line accross the screen
void PrintLine(char __c){ 
  for(int i = 0; i < 80; i++) cout << __c; 
  cout << endl;
}
//_____________________________________________________________________________

int main(int __argc,char *__argv[]){

  /*__________________________Parse the Command Line_________________________*/
  int c;
  // extern char* optarg;
  // extern int optind;
  
  while((c = getopt(__argc,__argv,"h")) != -1){
    switch(c){
    case 'h': // help option
      PrintUsage();
      return EXIT_SUCCESS;
      break;
    default:
      break;
    }
  }  
  /*____________________________Creating Matricies___________________________*/
  // Matrix is a template class, here we'll work with double's, but you 
  // could use any data type that defines the proper operators.
  Matrix<double> m22(2,2),m23(2,3),m33(3,3);
  IdentityMatrix<double> i2(2),i3(3);

  /*_____________________________Setting Matricies___________________________*/
  // The identity matrix is set defaultly when created. To set elements of our
  // other matricies we have a few options. We can either set them using 
  // other matricies (using the = operator), or by direct access to the 
  // elements.
  
  PrintLine(':');
  cout << "We've set up the following matricies: " << endl;
  // element access
  for(int i = 0; i < 2; i++) m22(i,i) = 2.0;
  m23(1,0) = -1.0;
  m23(0,2) = 2.0;
  m23(1,1) = 1.0;
  
  cout << "->m22:\n" << m22;
  cout << "->m23:\n" << m23;

  // using the assignment operator
  m33 = m23.T() * m23;
  cout << "->m33:\n" << m33;

  cout << "->identity(2x2):\n" << i2;
  cout << "->identity(3x3):\n" << i3;

  /*_____________________________Basic Operations____________________________*/
  PrintLine(':');
  cout << "Some basic matrix operations:" << endl;
  cout << "2*m22:\n" << 2*m22;
  cout << "m22 * m23:\n" << m22 * m23;
  cout << "m22 + m22:\n" << m22 + m22;
  cout << "trace(m33): " << m33.Trace() << endl;
  cout << "transpose(m23):\n" << m23.Transpose();
  cout << "etc...complete list of operations is in Matrix class documentation."
       << endl;

  PrintLine(':');

  return EXIT_SUCCESS;
}
//_____________________________________________________________________________

void PrintUsage(){ 
  cout << "Usage: matrix " << endl;
  cout << "This executable provides a number of example usages of the matrix "
       << "package. Run\nthe executable to see what's being done, then look at"
       << " the source file to see \nhow it's done in the code."
       << endl;
}
//_____________________________________________________________________________

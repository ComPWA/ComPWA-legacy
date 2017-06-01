#include "Examples/PythonDalitzFit/PythonFit.hpp"
#include <boost/python.hpp>

//using namespace ComPWA;

//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(X_returnsame_overloads, X::returnsame, 0, 1)
//BOOST_PYTHON_FUNCTION_OVERLOADS(X_returnsum_overloads, X::returnsum, 1, 2)
 
BOOST_PYTHON_MODULE(Dalitz_ext)
{
    using namespace boost::python;

    class_<PythonFit>("PythonFit")
    .def(init<PythonFit>())
    .def("StartFit", &PythonFit::StartFit)
    //.def("AddParameter", &ParameterList::AddParameter)
    //.def("to_str", &ParameterList::to_str)
    //.def("GetNDouble", &ParameterList::GetNDouble)
    ;

    //register_ptr_to_python< std::shared_ptr<Data> >();
    //register_ptr_to_python< std::shared_ptr<Amplitude> >();
    //register_ptr_to_python< std::shared_ptr<Optimizer::Optimizer> >();
   // register_ptr_to_python< std::shared_ptr<Generator> >();
}


/*USAGE
 * export COMPWA_DIR=YOUR_COMPWA_ROOT_DIR
 * export PYTHONPATH=$PYTHONPATH:YOUT_COMPWA_BUILD_DIR/lib
 * python
 * >>> import Dalitz_ext
 * >>> from Dalitz_ext import *
 * >>> fit = PythonFit()
 * >>> fit.StartFit()
 * >>> exit()
 */

#define PY_ARRAY_UNIQUE_SYMBOL CWRAPPER_ARRAY_API
#include "numpy_version.h"
#include "numpy/arrayobject.h"
#include "cwrapper.h"
#include <stdlib.h>
#include "ndarrays.h"
#include <stdint.h>
#include "cwrapper_ndarrays.h"


void bind_c_define_halosend(int64_t n0_w_c, double *w_c, int64_t n0_w_halosend, double *w_halosend, int64_t n0_indsend, int64_t *indsend);

/*........................................*/


/*........................................*/

/*........................................*/
PyObject *define_halosend_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray w_c = {.shape = NULL};
    t_ndarray w_halosend = {.shape = NULL};
    t_ndarray indsend = {.shape = NULL};
    PyArrayObject *w_c_tmp;
    PyArrayObject *w_halosend_tmp;
    PyArrayObject *indsend_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "w_c",
        "w_halosend",
        "indsend",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!", kwlist, &PyArray_Type, &w_c_tmp, &PyArray_Type, &w_halosend_tmp, &PyArray_Type, &indsend_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(w_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_c = pyarray_to_ndarray(w_c_tmp);
    }
    if (!pyarray_check(w_halosend_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_halosend = pyarray_to_ndarray(w_halosend_tmp);
    }
    if (!pyarray_check(indsend_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        indsend = pyarray_to_ndarray(indsend_tmp);
    }
    bind_c_define_halosend(nd_ndim(&w_c, 0), nd_data(&w_c), nd_ndim(&w_halosend, 0), nd_data(&w_halosend), nd_ndim(&indsend, 0), nd_data(&indsend));
    result = Py_BuildValue("");
    free_pointer(w_c);
    free_pointer(w_halosend);
    free_pointer(indsend);
    return result;
}
/*........................................*/

static int exec_func(PyObject* m)
{
    return 0;
}

/*........................................*/

static PyMethodDef pyccel_comm_methods[] = {
    {
        "define_halosend",
        (PyCFunction)define_halosend_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    { NULL, NULL, 0, NULL}
};

/*........................................*/

static PyModuleDef_Slot pyccel_comm_slots[] = {
    {Py_mod_exec, exec_func},
    {0, NULL},
};

/*........................................*/

static struct PyModuleDef pyccel_comm_module = {
    PyModuleDef_HEAD_INIT,
    /* name of module */
    "pyccel_comm",
    /* module documentation, may be NULL */
    NULL,
    /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    0,
    pyccel_comm_methods,
    pyccel_comm_slots
};

/*........................................*/

PyMODINIT_FUNC PyInit_pyccel_comm(void)
{
    import_array();
    return PyModuleDef_Init(&pyccel_comm_module);
}

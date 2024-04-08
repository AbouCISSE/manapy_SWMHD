#define PY_ARRAY_UNIQUE_SYMBOL CWRAPPER_ARRAY_API
#include "numpy_version.h"
#include "numpy/arrayobject.h"
#include "cwrapper.h"
#include <stdlib.h>
#include <stdint.h>
#include "ndarrays.h"
#include "cwrapper_ndarrays.h"


void bind_c_convert_solution(int64_t n0_x1, double *x1, int64_t n0_x1converted, double *x1converted, int64_t n0_tc, int64_t *tc, int64_t b0Size);
void bind_c_rhs_value_dirichlet_node(int64_t n0_Pbordnode, double *Pbordnode, int64_t n0_nodes, int64_t *nodes, int64_t n0_value, double *value);
void bind_c_rhs_value_dirichlet_face(int64_t n0_Pbordface, double *Pbordface, int64_t n0_faces, int64_t *faces, int64_t n0_value, double *value);
void bind_c_ghost_value_slip(int64_t n0_u_c, double *u_c, int64_t n0_v_c, double *v_c, int64_t n0_w_ghost, double *w_ghost, int64_t n0_cellid, int64_t n1_cellid, int64_t *cellid, int64_t n0_faces, int64_t *faces, int64_t n0_normal, int64_t n1_normal, double *normal, int64_t n0_mesure, double *mesure);
void bind_c_ghost_value_nonslip(int64_t n0_w_c, double *w_c, int64_t n0_w_ghost, double *w_ghost, int64_t n0_cellid, int64_t n1_cellid, int64_t *cellid, int64_t n0_faces, int64_t *faces);
void bind_c_ghost_value_neumann(int64_t n0_w_c, double *w_c, int64_t n0_w_ghost, double *w_ghost, int64_t n0_cellid, int64_t n1_cellid, int64_t *cellid, int64_t n0_faces, int64_t *faces);
void bind_c_ghost_value_dirichlet(int64_t n0_value, double *value, int64_t n0_w_ghost, double *w_ghost, int64_t n0_cellid, int64_t n1_cellid, int64_t *cellid, int64_t n0_faces, int64_t *faces);
void bind_c_haloghost_value_neumann(int64_t n0_w_halo, double *w_halo, int64_t n0_w_haloghost, double *w_haloghost, int64_t n0_haloghostcenter, int64_t n1_haloghostcenter, int64_t n2_haloghostcenter, double *haloghostcenter, int64_t BCindex, int64_t n0_halonodes, int64_t *halonodes);
void bind_c_haloghost_value_dirichlet(int64_t n0_value, double *value, int64_t n0_w_haloghost, double *w_haloghost, int64_t n0_haloghostcenter, int64_t n1_haloghostcenter, int64_t n2_haloghostcenter, double *haloghostcenter, int64_t BCindex, int64_t n0_halonodes, int64_t *halonodes);
void bind_c_haloghost_value_nonslip(int64_t n0_w_halo, double *w_halo, int64_t n0_w_haloghost, double *w_haloghost, int64_t n0_haloghostcenter, int64_t n1_haloghostcenter, int64_t n2_haloghostcenter, double *haloghostcenter, int64_t BCindex, int64_t n0_halonodes, int64_t *halonodes);
void bind_c_haloghost_value_slip(int64_t n0_u_halo, double *u_halo, int64_t n0_v_halo, double *v_halo, int64_t n0_w_haloghost, double *w_haloghost, int64_t n0_haloghostcenter, int64_t n1_haloghostcenter, int64_t n2_haloghostcenter, double *haloghostcenter, int64_t BCindex, int64_t n0_halonodes, int64_t *halonodes, int64_t n0_haloghostfaceinfo, int64_t n1_haloghostfaceinfo, int64_t n2_haloghostfaceinfo, double *haloghostfaceinfo);
void bind_c_cell_gradient_2d(int64_t n0_w_c, double *w_c, int64_t n0_w_ghost, double *w_ghost, int64_t n0_w_halo, double *w_halo, int64_t n0_w_haloghost, double *w_haloghost, int64_t n0_centerc, int64_t n1_centerc, double *centerc, int64_t n0_cellnid, int64_t n1_cellnid, int64_t *cellnid, int64_t n0_halonid, int64_t n1_halonid, int64_t *halonid, int64_t n0_nodecid, int64_t n1_nodecid, int64_t *nodecid, int64_t n0_periodicn, int64_t n1_periodicn, int64_t *periodicn, int64_t n0_periodic, int64_t n1_periodic, int64_t *periodic, int64_t n0_namen, int64_t *namen, int64_t n0_centerg, int64_t n1_centerg, int64_t n2_centerg, double *centerg, int64_t n0_halocenterg, int64_t n1_halocenterg, int64_t n2_halocenterg, double *halocenterg, int64_t n0_vertexn, int64_t n1_vertexn, double *vertexn, int64_t n0_centerh, int64_t n1_centerh, double *centerh, int64_t n0_shift, int64_t n1_shift, double *shift, int64_t nbproc, int64_t n0_w_x, double *w_x, int64_t n0_w_y, double *w_y, int64_t n0_w_z, double *w_z);
void bind_c_cell_gradient_3d(int64_t n0_w_c, double *w_c, int64_t n0_w_ghost, double *w_ghost, int64_t n0_w_halo, double *w_halo, int64_t n0_w_haloghost, double *w_haloghost, int64_t n0_centerc, int64_t n1_centerc, double *centerc, int64_t n0_cellnid, int64_t n1_cellnid, int64_t *cellnid, int64_t n0_halonid, int64_t n1_halonid, int64_t *halonid, int64_t n0_nodecid, int64_t n1_nodecid, int64_t *nodecid, int64_t n0_periodicn, int64_t n1_periodicn, int64_t *periodicn, int64_t n0_periodic, int64_t n1_periodic, int64_t *periodic, int64_t n0_namen, int64_t *namen, int64_t n0_centerg, int64_t n1_centerg, int64_t n2_centerg, double *centerg, int64_t n0_halocenterg, int64_t n1_halocenterg, int64_t n2_halocenterg, double *halocenterg, int64_t n0_vertexn, int64_t n1_vertexn, double *vertexn, int64_t n0_centerh, int64_t n1_centerh, double *centerh, int64_t n0_shift, int64_t n1_shift, double *shift, int64_t nbproc, int64_t n0_w_x, double *w_x, int64_t n0_w_y, double *w_y, int64_t n0_w_z, double *w_z);
void bind_c_face_gradient_2d(int64_t n0_w_c, double *w_c, int64_t n0_w_ghost, double *w_ghost, int64_t n0_w_halo, double *w_halo, int64_t n0_w_node, double *w_node, int64_t n0_cellidf, int64_t n1_cellidf, int64_t *cellidf, int64_t n0_nodeidf, int64_t n1_nodeidf, int64_t *nodeidf, int64_t n0_centergf, int64_t n1_centergf, double *centergf, int64_t n0_namef, int64_t *namef, int64_t n0_halofid, int64_t *halofid, int64_t n0_centerc, int64_t n1_centerc, double *centerc, int64_t n0_centerh, int64_t n1_centerh, double *centerh, int64_t n0_vertexn, int64_t n1_vertexn, double *vertexn, int64_t n0_airDiamond, double *airDiamond, int64_t n0_normalf, int64_t n1_normalf, double *normalf, int64_t n0_f_1, int64_t n1_f_1, double *f_1, int64_t n0_f_2, int64_t n1_f_2, double *f_2, int64_t n0_f_3, int64_t n1_f_3, double *f_3, int64_t n0_f_4, int64_t n1_f_4, double *f_4, int64_t n0_shift, int64_t n1_shift, double *shift, int64_t n0_wx_face, double *wx_face, int64_t n0_wy_face, double *wy_face, int64_t n0_wz_face, double *wz_face, int64_t n0_innerfaces, int64_t *innerfaces, int64_t n0_halofaces, int64_t *halofaces, int64_t n0_dirichletfaces, int64_t *dirichletfaces, int64_t n0_neumann, int64_t *neumann, int64_t n0_periodicfaces, int64_t *periodicfaces);
void bind_c_face_gradient_2d_uv(int64_t n0_h_c, double *h_c, int64_t n0_h_ghost, double *h_ghost, int64_t n0_h_halo, double *h_halo, int64_t n0_h_node, double *h_node, int64_t n0_w_c, double *w_c, int64_t n0_w_ghost, double *w_ghost, int64_t n0_w_halo, double *w_halo, int64_t n0_w_node, double *w_node, int64_t n0_cellidf, int64_t n1_cellidf, int64_t *cellidf, int64_t n0_nodeidf, int64_t n1_nodeidf, int64_t *nodeidf, int64_t n0_centergf, int64_t n1_centergf, double *centergf, int64_t n0_namef, int64_t *namef, int64_t n0_halofid, int64_t *halofid, int64_t n0_centerc, int64_t n1_centerc, double *centerc, int64_t n0_centerh, int64_t n1_centerh, double *centerh, int64_t n0_vertexn, int64_t n1_vertexn, double *vertexn, int64_t n0_airDiamond, double *airDiamond, int64_t n0_normalf, int64_t n1_normalf, double *normalf, int64_t n0_f_1, int64_t n1_f_1, double *f_1, int64_t n0_f_2, int64_t n1_f_2, double *f_2, int64_t n0_f_3, int64_t n1_f_3, double *f_3, int64_t n0_f_4, int64_t n1_f_4, double *f_4, int64_t n0_shift, int64_t n1_shift, double *shift, int64_t n0_wx_face, double *wx_face, int64_t n0_wy_face, double *wy_face, int64_t n0_wz_face, double *wz_face, int64_t n0_innerfaces, int64_t *innerfaces, int64_t n0_halofaces, int64_t *halofaces, int64_t n0_dirichletfaces, int64_t *dirichletfaces, int64_t n0_neumann, int64_t *neumann, int64_t n0_periodicfaces, int64_t *periodicfaces);
void bind_c_face_gradient_3d(int64_t n0_w_c, double *w_c, int64_t n0_w_ghost, double *w_ghost, int64_t n0_w_halo, double *w_halo, int64_t n0_w_node, double *w_node, int64_t n0_cellidf, int64_t n1_cellidf, int64_t *cellidf, int64_t n0_nodeidf, int64_t n1_nodeidf, int64_t *nodeidf, int64_t n0_centergf, int64_t n1_centergf, double *centergf, int64_t n0_namef, int64_t *namef, int64_t n0_halofid, int64_t *halofid, int64_t n0_centerc, int64_t n1_centerc, double *centerc, int64_t n0_centerh, int64_t n1_centerh, double *centerh, int64_t n0_vertexn, int64_t n1_vertexn, double *vertexn, int64_t n0_airDiamond, double *airDiamond, int64_t n0_normalf, int64_t n1_normalf, double *normalf, int64_t n0_f_1, int64_t n1_f_1, double *f_1, int64_t n0_f_2, int64_t n1_f_2, double *f_2, int64_t n0_f_3, int64_t n1_f_3, double *f_3, int64_t n0_f_4, int64_t n1_f_4, double *f_4, int64_t n0_shift, int64_t n1_shift, double *shift, int64_t n0_wx_face, double *wx_face, int64_t n0_wy_face, double *wy_face, int64_t n0_wz_face, double *wz_face, int64_t n0_innerfaces, int64_t *innerfaces, int64_t n0_halofaces, int64_t *halofaces, int64_t n0_dirichletfaces, int64_t *dirichletfaces, int64_t n0_neumann, int64_t *neumann, int64_t n0_periodicfaces, int64_t *periodicfaces);
void bind_c_centertovertex_2d(int64_t n0_w_c, double *w_c, int64_t n0_w_ghost, double *w_ghost, int64_t n0_w_halo, double *w_halo, int64_t n0_w_haloghost, double *w_haloghost, int64_t n0_centerc, int64_t n1_centerc, double *centerc, int64_t n0_centerh, int64_t n1_centerh, double *centerh, int64_t n0_cellidn, int64_t n1_cellidn, int64_t *cellidn, int64_t n0_periodicn, int64_t n1_periodicn, int64_t *periodicn, int64_t n0_haloidn, int64_t n1_haloidn, int64_t *haloidn, int64_t n0_vertexn, int64_t n1_vertexn, double *vertexn, int64_t n0_namen, int64_t *namen, int64_t n0_centergn, int64_t n1_centergn, int64_t n2_centergn, double *centergn, int64_t n0_halocentergn, int64_t n1_halocentergn, int64_t n2_halocentergn, double *halocentergn, int64_t n0_R_x, double *R_x, int64_t n0_R_y, double *R_y, int64_t n0_lambda_x, double *lambda_x, int64_t n0_lambda_y, double *lambda_y, int64_t n0_number, int64_t *number, int64_t n0_shift, int64_t n1_shift, double *shift, int64_t nbproc, int64_t n0_w_n, double *w_n);
void bind_c_centertovertex_3d(int64_t n0_w_c, double *w_c, int64_t n0_w_ghost, double *w_ghost, int64_t n0_w_halo, double *w_halo, int64_t n0_w_haloghost, double *w_haloghost, int64_t n0_centerc, int64_t n1_centerc, double *centerc, int64_t n0_centerh, int64_t n1_centerh, double *centerh, int64_t n0_cellidn, int64_t n1_cellidn, int64_t *cellidn, int64_t n0_periodicn, int64_t n1_periodicn, int64_t *periodicn, int64_t n0_haloidn, int64_t n1_haloidn, int64_t *haloidn, int64_t n0_vertexn, int64_t n1_vertexn, double *vertexn, int64_t n0_namen, int64_t *namen, int64_t n0_centergn, int64_t n1_centergn, int64_t n2_centergn, double *centergn, int64_t n0_halocentergn, int64_t n1_halocentergn, int64_t n2_halocentergn, double *halocentergn, int64_t n0_R_x, double *R_x, int64_t n0_R_y, double *R_y, int64_t n0_R_z, double *R_z, int64_t n0_lambda_x, double *lambda_x, int64_t n0_lambda_y, double *lambda_y, int64_t n0_lambda_z, double *lambda_z, int64_t n0_number, int64_t *number, int64_t n0_shift, int64_t n1_shift, double *shift, int64_t nbproc, int64_t n0_w_n, double *w_n);
void bind_c_barthlimiter_2d(int64_t n0_w_c, double *w_c, int64_t n0_w_ghost, double *w_ghost, int64_t n0_w_halo, double *w_halo, int64_t n0_w_x, double *w_x, int64_t n0_w_y, double *w_y, int64_t n0_w_z, double *w_z, int64_t n0_psi, double *psi, int64_t n0_cellid, int64_t n1_cellid, int64_t *cellid, int64_t n0_faceid, int64_t n1_faceid, int64_t *faceid, int64_t n0_namef, int64_t *namef, int64_t n0_halofid, int64_t *halofid, int64_t n0_centerc, int64_t n1_centerc, double *centerc, int64_t n0_centerf, int64_t n1_centerf, double *centerf);
void bind_c_barthlimiter_3d(int64_t n0_h_c, double *h_c, int64_t n0_h_ghost, double *h_ghost, int64_t n0_h_halo, double *h_halo, int64_t n0_h_x, double *h_x, int64_t n0_h_y, double *h_y, int64_t n0_h_z, double *h_z, int64_t n0_psi, double *psi, int64_t n0_cellid, int64_t n1_cellid, int64_t *cellid, int64_t n0_faceid, int64_t n1_faceid, int64_t *faceid, int64_t n0_namef, int64_t *namef, int64_t n0_halofid, int64_t *halofid, int64_t n0_centerc, int64_t n1_centerc, double *centerc, int64_t n0_centerf, int64_t n1_centerf, double *centerf);
int64_t bind_c_search_element(int64_t n0_a, int64_t *a, int64_t target_value);
void bind_c_get_triplet_2d(int64_t n0_cellfid, int64_t n1_cellfid, int64_t *cellfid, int64_t n0_nodeidf, int64_t n1_nodeidf, int64_t *nodeidf, int64_t n0_vertexn, int64_t n1_vertexn, double *vertexn, int64_t n0_halofid, int64_t *halofid, int64_t n0_haloext, int64_t n1_haloext, int64_t *haloext, int64_t n0_namen, int64_t *namen, int64_t n0_oldnamen, int64_t *oldnamen, int64_t n0_volume, double *volume, int64_t n0_cellnid, int64_t n1_cellnid, int64_t *cellnid, int64_t n0_centerc, int64_t n1_centerc, double *centerc, int64_t n0_centerh, int64_t n1_centerh, double *centerh, int64_t n0_halonid, int64_t n1_halonid, int64_t *halonid, int64_t n0_periodicnid, int64_t n1_periodicnid, int64_t *periodicnid, int64_t n0_centergn, int64_t n1_centergn, int64_t n2_centergn, double *centergn, int64_t n0_halocentergn, int64_t n1_halocentergn, int64_t n2_halocentergn, double *halocentergn, int64_t n0_airDiamond, double *airDiamond, int64_t n0_lambda_x, double *lambda_x, int64_t n0_lambda_y, double *lambda_y, int64_t n0_number, int64_t *number, int64_t n0_R_x, double *R_x, int64_t n0_R_y, double *R_y, int64_t n0_param1, double *param1, int64_t n0_param2, double *param2, int64_t n0_param3, double *param3, int64_t n0_param4, double *param4, int64_t n0_shift, int64_t n1_shift, double *shift, int64_t nbelements, int64_t n0_loctoglob, int64_t *loctoglob, int64_t n0_BCdirichlet, int64_t *BCdirichlet, int64_t n0_a_loc, double *a_loc, int64_t n0_irn_loc, int32_t *irn_loc, int64_t n0_jcn_loc, int32_t *jcn_loc, int64_t n0_matrixinnerfaces, int64_t *matrixinnerfaces, int64_t n0_halofaces, int64_t *halofaces, int64_t n0_dirichletfaces, int64_t *dirichletfaces);
int64_t bind_c_compute_2dmatrix_size(int64_t n0_nodeidf, int64_t n1_nodeidf, int64_t *nodeidf, int64_t n0_halofid, int64_t *halofid, int64_t n0_cellnid, int64_t n1_cellnid, int64_t *cellnid, int64_t n0_halonid, int64_t n1_halonid, int64_t *halonid, int64_t n0_periodicnid, int64_t n1_periodicnid, int64_t *periodicnid, int64_t n0_centergn, int64_t n1_centergn, int64_t n2_centergn, double *centergn, int64_t n0_halocentergn, int64_t n1_halocentergn, int64_t n2_halocentergn, double *halocentergn, int64_t n0_oldnamen, int64_t *oldnamen, int64_t n0_BCdirichlet, int64_t *BCdirichlet, int64_t n0_matrixinnerfaces, int64_t *matrixinnerfaces, int64_t n0_halofaces, int64_t *halofaces, int64_t n0_dirichletfaces, int64_t *dirichletfaces);
int64_t bind_c_compute_3dmatrix_size(int64_t n0_nodeidf, int64_t n1_nodeidf, int64_t *nodeidf, int64_t n0_halofid, int64_t *halofid, int64_t n0_cellnid, int64_t n1_cellnid, int64_t *cellnid, int64_t n0_halonid, int64_t n1_halonid, int64_t *halonid, int64_t n0_periodicnid, int64_t n1_periodicnid, int64_t *periodicnid, int64_t n0_centergn, int64_t n1_centergn, int64_t n2_centergn, double *centergn, int64_t n0_halocentergn, int64_t n1_halocentergn, int64_t n2_halocentergn, double *halocentergn, int64_t n0_oldnamen, int64_t *oldnamen, int64_t n0_BCdirichlet, int64_t *BCdirichlet, int64_t n0_matrixinnerfaces, int64_t *matrixinnerfaces, int64_t n0_halofaces, int64_t *halofaces, int64_t n0_dirichletfaces, int64_t *dirichletfaces);
void bind_c_get_triplet_3d(int64_t n0_cellfid, int64_t n1_cellfid, int64_t *cellfid, int64_t n0_nodeidf, int64_t n1_nodeidf, int64_t *nodeidf, int64_t n0_vertexn, int64_t n1_vertexn, double *vertexn, int64_t n0_halofid, int64_t *halofid, int64_t n0_haloext, int64_t n1_haloext, int64_t *haloext, int64_t n0_namen, int64_t *namen, int64_t n0_oldnamen, int64_t *oldnamen, int64_t n0_volume, double *volume, int64_t n0_centergn, int64_t n1_centergn, int64_t n2_centergn, double *centergn, int64_t n0_halocentergn, int64_t n1_halocentergn, int64_t n2_halocentergn, double *halocentergn, int64_t n0_periodicnid, int64_t n1_periodicnid, int64_t *periodicnid, int64_t n0_cellnid, int64_t n1_cellnid, int64_t *cellnid, int64_t n0_centerc, int64_t n1_centerc, double *centerc, int64_t n0_centerh, int64_t n1_centerh, double *centerh, int64_t n0_halonid, int64_t n1_halonid, int64_t *halonid, int64_t n0_airDiamond, double *airDiamond, int64_t n0_lambda_x, double *lambda_x, int64_t n0_lambda_y, double *lambda_y, int64_t n0_lambda_z, double *lambda_z, int64_t n0_number, int64_t *number, int64_t n0_R_x, double *R_x, int64_t n0_R_y, double *R_y, int64_t n0_R_z, double *R_z, int64_t n0_param1, double *param1, int64_t n0_param2, double *param2, int64_t n0_param3, double *param3, int64_t n0_shift, int64_t n1_shift, double *shift, int64_t n0_loctoglob, int64_t *loctoglob, int64_t n0_BCdirichlet, int64_t *BCdirichlet, int64_t n0_a_loc, double *a_loc, int64_t n0_irn_loc, int32_t *irn_loc, int64_t n0_jcn_loc, int32_t *jcn_loc, int64_t n0_matrixinnerfaces, int64_t *matrixinnerfaces, int64_t n0_halofaces, int64_t *halofaces, int64_t n0_dirichletfaces, int64_t *dirichletfaces);
void bind_c_get_rhs_loc_2d(int64_t n0_cellfid, int64_t n1_cellfid, int64_t *cellfid, int64_t n0_nodeidf, int64_t n1_nodeidf, int64_t *nodeidf, int64_t n0_oldname, int64_t *oldname, int64_t n0_volume, double *volume, int64_t n0_centergn, int64_t n1_centergn, int64_t n2_centergn, double *centergn, int64_t n0_param1, double *param1, int64_t n0_param2, double *param2, int64_t n0_param3, double *param3, int64_t n0_param4, double *param4, int64_t n0_Pbordnode, double *Pbordnode, int64_t n0_Pbordface, double *Pbordface, int64_t n0_rhs_loc, double *rhs_loc, int64_t n0_BCdirichlet, int64_t *BCdirichlet, int64_t n0_centergf, int64_t n1_centergf, double *centergf, int64_t n0_matrixinnerfaces, int64_t *matrixinnerfaces, int64_t n0_halofaces, int64_t *halofaces, int64_t n0_dirichletfaces, int64_t *dirichletfaces);
void bind_c_get_rhs_glob_2d(int64_t n0_cellfid, int64_t n1_cellfid, int64_t *cellfid, int64_t n0_nodeidf, int64_t n1_nodeidf, int64_t *nodeidf, int64_t n0_oldname, int64_t *oldname, int64_t n0_volume, double *volume, int64_t n0_centergn, int64_t n1_centergn, int64_t n2_centergn, double *centergn, int64_t n0_loctoglob, int64_t *loctoglob, int64_t n0_param1, double *param1, int64_t n0_param2, double *param2, int64_t n0_param3, double *param3, int64_t n0_param4, double *param4, int64_t n0_Pbordnode, double *Pbordnode, int64_t n0_Pbordface, double *Pbordface, int64_t n0_rhs, double *rhs, int64_t n0_BCdirichlet, int64_t *BCdirichlet, int64_t n0_centergf, int64_t n1_centergf, double *centergf, int64_t n0_matrixinnerfaces, int64_t *matrixinnerfaces, int64_t n0_halofaces, int64_t *halofaces, int64_t n0_dirichletfaces, int64_t *dirichletfaces);
void bind_c_get_rhs_loc_3d(int64_t n0_cellfid, int64_t n1_cellfid, int64_t *cellfid, int64_t n0_nodeidf, int64_t n1_nodeidf, int64_t *nodeidf, int64_t n0_oldname, int64_t *oldname, int64_t n0_volume, double *volume, int64_t n0_centergn, int64_t n1_centergn, int64_t n2_centergn, double *centergn, int64_t n0_param1, double *param1, int64_t n0_param2, double *param2, int64_t n0_param3, double *param3, int64_t n0_Pbordnode, double *Pbordnode, int64_t n0_Pbordface, double *Pbordface, int64_t n0_rhs_loc, double *rhs_loc, int64_t n0_BCdirichlet, int64_t *BCdirichlet, int64_t n0_matrixinnerfaces, int64_t *matrixinnerfaces, int64_t n0_halofaces, int64_t *halofaces, int64_t n0_dirichletfaces, int64_t *dirichletfaces);
void bind_c_get_rhs_glob_3d(int64_t n0_cellfid, int64_t n1_cellfid, int64_t *cellfid, int64_t n0_nodeidf, int64_t n1_nodeidf, int64_t *nodeidf, int64_t n0_oldname, int64_t *oldname, int64_t n0_volume, double *volume, int64_t n0_centergn, int64_t n1_centergn, int64_t n2_centergn, double *centergn, int64_t n0_loctoglob, int64_t *loctoglob, int64_t n0_param1, double *param1, int64_t n0_param2, double *param2, int64_t n0_param3, double *param3, int64_t n0_Pbordnode, double *Pbordnode, int64_t n0_Pbordface, double *Pbordface, int64_t n0_rhs, double *rhs, int64_t n0_BCdirichlet, int64_t *BCdirichlet, int64_t n0_matrixinnerfaces, int64_t *matrixinnerfaces, int64_t n0_halofaces, int64_t *halofaces, int64_t n0_dirichletfaces, int64_t *dirichletfaces);
void bind_c_compute_p_gradient_2d(int64_t n0_P_c, double *P_c, int64_t n0_P_ghost, double *P_ghost, int64_t n0_P_halo, double *P_halo, int64_t n0_P_node, double *P_node, int64_t n0_cellidf, int64_t n1_cellidf, int64_t *cellidf, int64_t n0_nodeidf, int64_t n1_nodeidf, int64_t *nodeidf, int64_t n0_centergf, int64_t n1_centergf, double *centergf, int64_t n0_namef, int64_t *namef, int64_t n0_halofid, int64_t *halofid, int64_t n0_centerc, int64_t n1_centerc, double *centerc, int64_t n0_centerh, int64_t n1_centerh, double *centerh, int64_t n0_oldname, int64_t *oldname, int64_t n0_airDiamond, double *airDiamond, int64_t n0_f_1, int64_t n1_f_1, double *f_1, int64_t n0_f_2, int64_t n1_f_2, double *f_2, int64_t n0_f_3, int64_t n1_f_3, double *f_3, int64_t n0_f_4, int64_t n1_f_4, double *f_4, int64_t n0_normalf, int64_t n1_normalf, double *normalf, int64_t n0_shift, int64_t n1_shift, double *shift, int64_t n0_Pbordnode, double *Pbordnode, int64_t n0_Pbordface, double *Pbordface, int64_t n0_Px_face, double *Px_face, int64_t n0_Py_face, double *Py_face, int64_t n0_Pz_face, double *Pz_face, int64_t n0_BCdirichlet, int64_t *BCdirichlet, int64_t n0_innerfaces, int64_t *innerfaces, int64_t n0_halofaces, int64_t *halofaces, int64_t n0_neumannfaces, int64_t *neumannfaces, int64_t n0_dirichletfaces, int64_t *dirichletfaces, int64_t n0_periodicfaces, int64_t *periodicfaces);
void bind_c_compute_p_gradient_3d(int64_t n0_val_c, double *val_c, int64_t n0_v_ghost, double *v_ghost, int64_t n0_v_halo, double *v_halo, int64_t n0_v_node, double *v_node, int64_t n0_cellidf, int64_t n1_cellidf, int64_t *cellidf, int64_t n0_nodeidf, int64_t n1_nodeidf, int64_t *nodeidf, int64_t n0_centergf, int64_t n1_centergf, double *centergf, int64_t n0_namef, int64_t *namef, int64_t n0_halofid, int64_t *halofid, int64_t n0_centerc, int64_t n1_centerc, double *centerc, int64_t n0_centerh, int64_t n1_centerh, double *centerh, int64_t n0_oldname, int64_t *oldname, int64_t n0_airDiamond, double *airDiamond, int64_t n0_n1, int64_t n1_n1, double *n1, int64_t n0_n2, int64_t n1_n2, double *n2, int64_t n0_n3, int64_t n1_n3, double *n3, int64_t n0_n4, int64_t n1_n4, double *n4, int64_t n0_normalf, int64_t n1_normalf, double *normalf, int64_t n0_shift, int64_t n1_shift, double *shift, int64_t n0_Pbordnode, double *Pbordnode, int64_t n0_Pbordface, double *Pbordface, int64_t n0_Px_face, double *Px_face, int64_t n0_Py_face, double *Py_face, int64_t n0_Pz_face, double *Pz_face, int64_t n0_BCdirichlet, int64_t *BCdirichlet, int64_t n0_innerfaces, int64_t *innerfaces, int64_t n0_halofaces, int64_t *halofaces, int64_t n0_neumannfaces, int64_t *neumannfaces, int64_t n0_dirichletfaces, int64_t *dirichletfaces, int64_t n0_periodicfaces, int64_t *periodicfaces);
void bind_c_facetocell(int64_t n0_u_face, double *u_face, int64_t n0_u_c, double *u_c, int64_t n0_faceidc, int64_t n1_faceidc, int64_t *faceidc, int64_t dim);

/*........................................*/


/*........................................*/

/*........................................*/
PyObject *convert_solution_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray x1 = {.shape = NULL};
    t_ndarray x1converted = {.shape = NULL};
    t_ndarray tc = {.shape = NULL};
    int64_t b0Size;
    PyArrayObject *x1_tmp;
    PyArrayObject *x1converted_tmp;
    PyArrayObject *tc_tmp;
    PyObject *b0Size_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "x1",
        "x1converted",
        "tc",
        "b0Size",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O", kwlist, &PyArray_Type, &x1_tmp, &PyArray_Type, &x1converted_tmp, &PyArray_Type, &tc_tmp, &b0Size_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(x1_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        x1 = pyarray_to_ndarray(x1_tmp);
    }
    if (!pyarray_check(x1converted_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        x1converted = pyarray_to_ndarray(x1converted_tmp);
    }
    if (!pyarray_check(tc_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        tc = pyarray_to_ndarray(tc_tmp);
    }
    if (PyIs_NativeInt(b0Size_tmp))
    {
        b0Size = PyInt64_to_Int64(b0Size_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    bind_c_convert_solution(nd_ndim(&x1, 0), nd_data(&x1), nd_ndim(&x1converted, 0), nd_data(&x1converted), nd_ndim(&tc, 0), nd_data(&tc), b0Size);
    result = Py_BuildValue("");
    free_pointer(x1);
    free_pointer(x1converted);
    free_pointer(tc);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *rhs_value_dirichlet_node_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray Pbordnode = {.shape = NULL};
    t_ndarray nodes = {.shape = NULL};
    t_ndarray value = {.shape = NULL};
    PyArrayObject *Pbordnode_tmp;
    PyArrayObject *nodes_tmp;
    PyArrayObject *value_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "Pbordnode",
        "nodes",
        "value",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!", kwlist, &PyArray_Type, &Pbordnode_tmp, &PyArray_Type, &nodes_tmp, &PyArray_Type, &value_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(Pbordnode_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        Pbordnode = pyarray_to_ndarray(Pbordnode_tmp);
    }
    if (!pyarray_check(nodes_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        nodes = pyarray_to_ndarray(nodes_tmp);
    }
    if (!pyarray_check(value_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        value = pyarray_to_ndarray(value_tmp);
    }
    bind_c_rhs_value_dirichlet_node(nd_ndim(&Pbordnode, 0), nd_data(&Pbordnode), nd_ndim(&nodes, 0), nd_data(&nodes), nd_ndim(&value, 0), nd_data(&value));
    result = Py_BuildValue("");
    free_pointer(Pbordnode);
    free_pointer(nodes);
    free_pointer(value);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *rhs_value_dirichlet_face_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray Pbordface = {.shape = NULL};
    t_ndarray faces = {.shape = NULL};
    t_ndarray value = {.shape = NULL};
    PyArrayObject *Pbordface_tmp;
    PyArrayObject *faces_tmp;
    PyArrayObject *value_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "Pbordface",
        "faces",
        "value",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!", kwlist, &PyArray_Type, &Pbordface_tmp, &PyArray_Type, &faces_tmp, &PyArray_Type, &value_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(Pbordface_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        Pbordface = pyarray_to_ndarray(Pbordface_tmp);
    }
    if (!pyarray_check(faces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        faces = pyarray_to_ndarray(faces_tmp);
    }
    if (!pyarray_check(value_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        value = pyarray_to_ndarray(value_tmp);
    }
    bind_c_rhs_value_dirichlet_face(nd_ndim(&Pbordface, 0), nd_data(&Pbordface), nd_ndim(&faces, 0), nd_data(&faces), nd_ndim(&value, 0), nd_data(&value));
    result = Py_BuildValue("");
    free_pointer(Pbordface);
    free_pointer(faces);
    free_pointer(value);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *ghost_value_slip_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray u_c = {.shape = NULL};
    t_ndarray v_c = {.shape = NULL};
    t_ndarray w_ghost = {.shape = NULL};
    t_ndarray cellid = {.shape = NULL};
    t_ndarray faces = {.shape = NULL};
    t_ndarray normal = {.shape = NULL};
    t_ndarray mesure = {.shape = NULL};
    PyArrayObject *u_c_tmp;
    PyArrayObject *v_c_tmp;
    PyArrayObject *w_ghost_tmp;
    PyArrayObject *cellid_tmp;
    PyArrayObject *faces_tmp;
    PyArrayObject *normal_tmp;
    PyArrayObject *mesure_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "u_c",
        "v_c",
        "w_ghost",
        "cellid",
        "faces",
        "normal",
        "mesure",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!", kwlist, &PyArray_Type, &u_c_tmp, &PyArray_Type, &v_c_tmp, &PyArray_Type, &w_ghost_tmp, &PyArray_Type, &cellid_tmp, &PyArray_Type, &faces_tmp, &PyArray_Type, &normal_tmp, &PyArray_Type, &mesure_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(u_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        u_c = pyarray_to_ndarray(u_c_tmp);
    }
    if (!pyarray_check(v_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        v_c = pyarray_to_ndarray(v_c_tmp);
    }
    if (!pyarray_check(w_ghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_ghost = pyarray_to_ndarray(w_ghost_tmp);
    }
    if (!pyarray_check(cellid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellid = pyarray_to_ndarray(cellid_tmp);
    }
    if (!pyarray_check(faces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        faces = pyarray_to_ndarray(faces_tmp);
    }
    if (!pyarray_check(normal_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        normal = pyarray_to_ndarray(normal_tmp);
    }
    if (!pyarray_check(mesure_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        mesure = pyarray_to_ndarray(mesure_tmp);
    }
    bind_c_ghost_value_slip(nd_ndim(&u_c, 0), nd_data(&u_c), nd_ndim(&v_c, 0), nd_data(&v_c), nd_ndim(&w_ghost, 0), nd_data(&w_ghost), nd_ndim(&cellid, 0), nd_ndim(&cellid, 1), nd_data(&cellid), nd_ndim(&faces, 0), nd_data(&faces), nd_ndim(&normal, 0), nd_ndim(&normal, 1), nd_data(&normal), nd_ndim(&mesure, 0), nd_data(&mesure));
    result = Py_BuildValue("");
    free_pointer(u_c);
    free_pointer(v_c);
    free_pointer(w_ghost);
    free_pointer(cellid);
    free_pointer(faces);
    free_pointer(normal);
    free_pointer(mesure);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *ghost_value_nonslip_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray w_c = {.shape = NULL};
    t_ndarray w_ghost = {.shape = NULL};
    t_ndarray cellid = {.shape = NULL};
    t_ndarray faces = {.shape = NULL};
    PyArrayObject *w_c_tmp;
    PyArrayObject *w_ghost_tmp;
    PyArrayObject *cellid_tmp;
    PyArrayObject *faces_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "w_c",
        "w_ghost",
        "cellid",
        "faces",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!", kwlist, &PyArray_Type, &w_c_tmp, &PyArray_Type, &w_ghost_tmp, &PyArray_Type, &cellid_tmp, &PyArray_Type, &faces_tmp))
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
    if (!pyarray_check(w_ghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_ghost = pyarray_to_ndarray(w_ghost_tmp);
    }
    if (!pyarray_check(cellid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellid = pyarray_to_ndarray(cellid_tmp);
    }
    if (!pyarray_check(faces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        faces = pyarray_to_ndarray(faces_tmp);
    }
    bind_c_ghost_value_nonslip(nd_ndim(&w_c, 0), nd_data(&w_c), nd_ndim(&w_ghost, 0), nd_data(&w_ghost), nd_ndim(&cellid, 0), nd_ndim(&cellid, 1), nd_data(&cellid), nd_ndim(&faces, 0), nd_data(&faces));
    result = Py_BuildValue("");
    free_pointer(w_c);
    free_pointer(w_ghost);
    free_pointer(cellid);
    free_pointer(faces);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *ghost_value_neumann_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray w_c = {.shape = NULL};
    t_ndarray w_ghost = {.shape = NULL};
    t_ndarray cellid = {.shape = NULL};
    t_ndarray faces = {.shape = NULL};
    PyArrayObject *w_c_tmp;
    PyArrayObject *w_ghost_tmp;
    PyArrayObject *cellid_tmp;
    PyArrayObject *faces_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "w_c",
        "w_ghost",
        "cellid",
        "faces",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!", kwlist, &PyArray_Type, &w_c_tmp, &PyArray_Type, &w_ghost_tmp, &PyArray_Type, &cellid_tmp, &PyArray_Type, &faces_tmp))
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
    if (!pyarray_check(w_ghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_ghost = pyarray_to_ndarray(w_ghost_tmp);
    }
    if (!pyarray_check(cellid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellid = pyarray_to_ndarray(cellid_tmp);
    }
    if (!pyarray_check(faces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        faces = pyarray_to_ndarray(faces_tmp);
    }
    bind_c_ghost_value_neumann(nd_ndim(&w_c, 0), nd_data(&w_c), nd_ndim(&w_ghost, 0), nd_data(&w_ghost), nd_ndim(&cellid, 0), nd_ndim(&cellid, 1), nd_data(&cellid), nd_ndim(&faces, 0), nd_data(&faces));
    result = Py_BuildValue("");
    free_pointer(w_c);
    free_pointer(w_ghost);
    free_pointer(cellid);
    free_pointer(faces);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *ghost_value_dirichlet_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray value = {.shape = NULL};
    t_ndarray w_ghost = {.shape = NULL};
    t_ndarray cellid = {.shape = NULL};
    t_ndarray faces = {.shape = NULL};
    PyArrayObject *value_tmp;
    PyArrayObject *w_ghost_tmp;
    PyArrayObject *cellid_tmp;
    PyArrayObject *faces_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "value",
        "w_ghost",
        "cellid",
        "faces",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!", kwlist, &PyArray_Type, &value_tmp, &PyArray_Type, &w_ghost_tmp, &PyArray_Type, &cellid_tmp, &PyArray_Type, &faces_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(value_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        value = pyarray_to_ndarray(value_tmp);
    }
    if (!pyarray_check(w_ghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_ghost = pyarray_to_ndarray(w_ghost_tmp);
    }
    if (!pyarray_check(cellid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellid = pyarray_to_ndarray(cellid_tmp);
    }
    if (!pyarray_check(faces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        faces = pyarray_to_ndarray(faces_tmp);
    }
    bind_c_ghost_value_dirichlet(nd_ndim(&value, 0), nd_data(&value), nd_ndim(&w_ghost, 0), nd_data(&w_ghost), nd_ndim(&cellid, 0), nd_ndim(&cellid, 1), nd_data(&cellid), nd_ndim(&faces, 0), nd_data(&faces));
    result = Py_BuildValue("");
    free_pointer(value);
    free_pointer(w_ghost);
    free_pointer(cellid);
    free_pointer(faces);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *haloghost_value_neumann_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray w_halo = {.shape = NULL};
    t_ndarray w_haloghost = {.shape = NULL};
    t_ndarray haloghostcenter = {.shape = NULL};
    int64_t BCindex;
    t_ndarray halonodes = {.shape = NULL};
    PyArrayObject *w_halo_tmp;
    PyArrayObject *w_haloghost_tmp;
    PyArrayObject *haloghostcenter_tmp;
    PyObject *BCindex_tmp;
    PyArrayObject *halonodes_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "w_halo",
        "w_haloghost",
        "haloghostcenter",
        "BCindex",
        "halonodes",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!OO!", kwlist, &PyArray_Type, &w_halo_tmp, &PyArray_Type, &w_haloghost_tmp, &PyArray_Type, &haloghostcenter_tmp, &BCindex_tmp, &PyArray_Type, &halonodes_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(w_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_halo = pyarray_to_ndarray(w_halo_tmp);
    }
    if (!pyarray_check(w_haloghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_haloghost = pyarray_to_ndarray(w_haloghost_tmp);
    }
    if (!pyarray_check(haloghostcenter_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        haloghostcenter = pyarray_to_ndarray(haloghostcenter_tmp);
    }
    if (PyIs_NativeInt(BCindex_tmp))
    {
        BCindex = PyInt64_to_Int64(BCindex_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    if (!pyarray_check(halonodes_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halonodes = pyarray_to_ndarray(halonodes_tmp);
    }
    bind_c_haloghost_value_neumann(nd_ndim(&w_halo, 0), nd_data(&w_halo), nd_ndim(&w_haloghost, 0), nd_data(&w_haloghost), nd_ndim(&haloghostcenter, 0), nd_ndim(&haloghostcenter, 1), nd_ndim(&haloghostcenter, 2), nd_data(&haloghostcenter), BCindex, nd_ndim(&halonodes, 0), nd_data(&halonodes));
    result = Py_BuildValue("");
    free_pointer(w_halo);
    free_pointer(w_haloghost);
    free_pointer(haloghostcenter);
    free_pointer(halonodes);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *haloghost_value_dirichlet_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray value = {.shape = NULL};
    t_ndarray w_haloghost = {.shape = NULL};
    t_ndarray haloghostcenter = {.shape = NULL};
    int64_t BCindex;
    t_ndarray halonodes = {.shape = NULL};
    PyArrayObject *value_tmp;
    PyArrayObject *w_haloghost_tmp;
    PyArrayObject *haloghostcenter_tmp;
    PyObject *BCindex_tmp;
    PyArrayObject *halonodes_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "value",
        "w_haloghost",
        "haloghostcenter",
        "BCindex",
        "halonodes",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!OO!", kwlist, &PyArray_Type, &value_tmp, &PyArray_Type, &w_haloghost_tmp, &PyArray_Type, &haloghostcenter_tmp, &BCindex_tmp, &PyArray_Type, &halonodes_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(value_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        value = pyarray_to_ndarray(value_tmp);
    }
    if (!pyarray_check(w_haloghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_haloghost = pyarray_to_ndarray(w_haloghost_tmp);
    }
    if (!pyarray_check(haloghostcenter_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        haloghostcenter = pyarray_to_ndarray(haloghostcenter_tmp);
    }
    if (PyIs_NativeInt(BCindex_tmp))
    {
        BCindex = PyInt64_to_Int64(BCindex_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    if (!pyarray_check(halonodes_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halonodes = pyarray_to_ndarray(halonodes_tmp);
    }
    bind_c_haloghost_value_dirichlet(nd_ndim(&value, 0), nd_data(&value), nd_ndim(&w_haloghost, 0), nd_data(&w_haloghost), nd_ndim(&haloghostcenter, 0), nd_ndim(&haloghostcenter, 1), nd_ndim(&haloghostcenter, 2), nd_data(&haloghostcenter), BCindex, nd_ndim(&halonodes, 0), nd_data(&halonodes));
    result = Py_BuildValue("");
    free_pointer(value);
    free_pointer(w_haloghost);
    free_pointer(haloghostcenter);
    free_pointer(halonodes);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *haloghost_value_nonslip_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray w_halo = {.shape = NULL};
    t_ndarray w_haloghost = {.shape = NULL};
    t_ndarray haloghostcenter = {.shape = NULL};
    int64_t BCindex;
    t_ndarray halonodes = {.shape = NULL};
    PyArrayObject *w_halo_tmp;
    PyArrayObject *w_haloghost_tmp;
    PyArrayObject *haloghostcenter_tmp;
    PyObject *BCindex_tmp;
    PyArrayObject *halonodes_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "w_halo",
        "w_haloghost",
        "haloghostcenter",
        "BCindex",
        "halonodes",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!OO!", kwlist, &PyArray_Type, &w_halo_tmp, &PyArray_Type, &w_haloghost_tmp, &PyArray_Type, &haloghostcenter_tmp, &BCindex_tmp, &PyArray_Type, &halonodes_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(w_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_halo = pyarray_to_ndarray(w_halo_tmp);
    }
    if (!pyarray_check(w_haloghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_haloghost = pyarray_to_ndarray(w_haloghost_tmp);
    }
    if (!pyarray_check(haloghostcenter_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        haloghostcenter = pyarray_to_ndarray(haloghostcenter_tmp);
    }
    if (PyIs_NativeInt(BCindex_tmp))
    {
        BCindex = PyInt64_to_Int64(BCindex_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    if (!pyarray_check(halonodes_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halonodes = pyarray_to_ndarray(halonodes_tmp);
    }
    bind_c_haloghost_value_nonslip(nd_ndim(&w_halo, 0), nd_data(&w_halo), nd_ndim(&w_haloghost, 0), nd_data(&w_haloghost), nd_ndim(&haloghostcenter, 0), nd_ndim(&haloghostcenter, 1), nd_ndim(&haloghostcenter, 2), nd_data(&haloghostcenter), BCindex, nd_ndim(&halonodes, 0), nd_data(&halonodes));
    result = Py_BuildValue("");
    free_pointer(w_halo);
    free_pointer(w_haloghost);
    free_pointer(haloghostcenter);
    free_pointer(halonodes);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *haloghost_value_slip_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray u_halo = {.shape = NULL};
    t_ndarray v_halo = {.shape = NULL};
    t_ndarray w_haloghost = {.shape = NULL};
    t_ndarray haloghostcenter = {.shape = NULL};
    int64_t BCindex;
    t_ndarray halonodes = {.shape = NULL};
    t_ndarray haloghostfaceinfo = {.shape = NULL};
    PyArrayObject *u_halo_tmp;
    PyArrayObject *v_halo_tmp;
    PyArrayObject *w_haloghost_tmp;
    PyArrayObject *haloghostcenter_tmp;
    PyObject *BCindex_tmp;
    PyArrayObject *halonodes_tmp;
    PyArrayObject *haloghostfaceinfo_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "u_halo",
        "v_halo",
        "w_haloghost",
        "haloghostcenter",
        "BCindex",
        "halonodes",
        "haloghostfaceinfo",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!OO!O!", kwlist, &PyArray_Type, &u_halo_tmp, &PyArray_Type, &v_halo_tmp, &PyArray_Type, &w_haloghost_tmp, &PyArray_Type, &haloghostcenter_tmp, &BCindex_tmp, &PyArray_Type, &halonodes_tmp, &PyArray_Type, &haloghostfaceinfo_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(u_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        u_halo = pyarray_to_ndarray(u_halo_tmp);
    }
    if (!pyarray_check(v_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        v_halo = pyarray_to_ndarray(v_halo_tmp);
    }
    if (!pyarray_check(w_haloghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_haloghost = pyarray_to_ndarray(w_haloghost_tmp);
    }
    if (!pyarray_check(haloghostcenter_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        haloghostcenter = pyarray_to_ndarray(haloghostcenter_tmp);
    }
    if (PyIs_NativeInt(BCindex_tmp))
    {
        BCindex = PyInt64_to_Int64(BCindex_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    if (!pyarray_check(halonodes_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halonodes = pyarray_to_ndarray(halonodes_tmp);
    }
    if (!pyarray_check(haloghostfaceinfo_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        haloghostfaceinfo = pyarray_to_ndarray(haloghostfaceinfo_tmp);
    }
    bind_c_haloghost_value_slip(nd_ndim(&u_halo, 0), nd_data(&u_halo), nd_ndim(&v_halo, 0), nd_data(&v_halo), nd_ndim(&w_haloghost, 0), nd_data(&w_haloghost), nd_ndim(&haloghostcenter, 0), nd_ndim(&haloghostcenter, 1), nd_ndim(&haloghostcenter, 2), nd_data(&haloghostcenter), BCindex, nd_ndim(&halonodes, 0), nd_data(&halonodes), nd_ndim(&haloghostfaceinfo, 0), nd_ndim(&haloghostfaceinfo, 1), nd_ndim(&haloghostfaceinfo, 2), nd_data(&haloghostfaceinfo));
    result = Py_BuildValue("");
    free_pointer(u_halo);
    free_pointer(v_halo);
    free_pointer(w_haloghost);
    free_pointer(haloghostcenter);
    free_pointer(halonodes);
    free_pointer(haloghostfaceinfo);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *cell_gradient_2d_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray w_c = {.shape = NULL};
    t_ndarray w_ghost = {.shape = NULL};
    t_ndarray w_halo = {.shape = NULL};
    t_ndarray w_haloghost = {.shape = NULL};
    t_ndarray centerc = {.shape = NULL};
    t_ndarray cellnid = {.shape = NULL};
    t_ndarray halonid = {.shape = NULL};
    t_ndarray nodecid = {.shape = NULL};
    t_ndarray periodicn = {.shape = NULL};
    t_ndarray periodic = {.shape = NULL};
    t_ndarray namen = {.shape = NULL};
    t_ndarray centerg = {.shape = NULL};
    t_ndarray halocenterg = {.shape = NULL};
    t_ndarray vertexn = {.shape = NULL};
    t_ndarray centerh = {.shape = NULL};
    t_ndarray shift = {.shape = NULL};
    int64_t nbproc;
    t_ndarray w_x = {.shape = NULL};
    t_ndarray w_y = {.shape = NULL};
    t_ndarray w_z = {.shape = NULL};
    PyArrayObject *w_c_tmp;
    PyArrayObject *w_ghost_tmp;
    PyArrayObject *w_halo_tmp;
    PyArrayObject *w_haloghost_tmp;
    PyArrayObject *centerc_tmp;
    PyArrayObject *cellnid_tmp;
    PyArrayObject *halonid_tmp;
    PyArrayObject *nodecid_tmp;
    PyArrayObject *periodicn_tmp;
    PyArrayObject *periodic_tmp;
    PyArrayObject *namen_tmp;
    PyArrayObject *centerg_tmp;
    PyArrayObject *halocenterg_tmp;
    PyArrayObject *vertexn_tmp;
    PyArrayObject *centerh_tmp;
    PyArrayObject *shift_tmp;
    PyObject *nbproc_tmp;
    PyArrayObject *w_x_tmp;
    PyArrayObject *w_y_tmp;
    PyArrayObject *w_z_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "w_c",
        "w_ghost",
        "w_halo",
        "w_haloghost",
        "centerc",
        "cellnid",
        "halonid",
        "nodecid",
        "periodicn",
        "periodic",
        "namen",
        "centerg",
        "halocenterg",
        "vertexn",
        "centerh",
        "shift",
        "nbproc",
        "w_x",
        "w_y",
        "w_z",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!OO!O!O!", kwlist, &PyArray_Type, &w_c_tmp, &PyArray_Type, &w_ghost_tmp, &PyArray_Type, &w_halo_tmp, &PyArray_Type, &w_haloghost_tmp, &PyArray_Type, &centerc_tmp, &PyArray_Type, &cellnid_tmp, &PyArray_Type, &halonid_tmp, &PyArray_Type, &nodecid_tmp, &PyArray_Type, &periodicn_tmp, &PyArray_Type, &periodic_tmp, &PyArray_Type, &namen_tmp, &PyArray_Type, &centerg_tmp, &PyArray_Type, &halocenterg_tmp, &PyArray_Type, &vertexn_tmp, &PyArray_Type, &centerh_tmp, &PyArray_Type, &shift_tmp, &nbproc_tmp, &PyArray_Type, &w_x_tmp, &PyArray_Type, &w_y_tmp, &PyArray_Type, &w_z_tmp))
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
    if (!pyarray_check(w_ghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_ghost = pyarray_to_ndarray(w_ghost_tmp);
    }
    if (!pyarray_check(w_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_halo = pyarray_to_ndarray(w_halo_tmp);
    }
    if (!pyarray_check(w_haloghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_haloghost = pyarray_to_ndarray(w_haloghost_tmp);
    }
    if (!pyarray_check(centerc_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerc = pyarray_to_ndarray(centerc_tmp);
    }
    if (!pyarray_check(cellnid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellnid = pyarray_to_ndarray(cellnid_tmp);
    }
    if (!pyarray_check(halonid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        halonid = pyarray_to_ndarray(halonid_tmp);
    }
    if (!pyarray_check(nodecid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nodecid = pyarray_to_ndarray(nodecid_tmp);
    }
    if (!pyarray_check(periodicn_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        periodicn = pyarray_to_ndarray(periodicn_tmp);
    }
    if (!pyarray_check(periodic_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        periodic = pyarray_to_ndarray(periodic_tmp);
    }
    if (!pyarray_check(namen_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        namen = pyarray_to_ndarray(namen_tmp);
    }
    if (!pyarray_check(centerg_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerg = pyarray_to_ndarray(centerg_tmp);
    }
    if (!pyarray_check(halocenterg_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        halocenterg = pyarray_to_ndarray(halocenterg_tmp);
    }
    if (!pyarray_check(vertexn_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        vertexn = pyarray_to_ndarray(vertexn_tmp);
    }
    if (!pyarray_check(centerh_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerh = pyarray_to_ndarray(centerh_tmp);
    }
    if (!pyarray_check(shift_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        shift = pyarray_to_ndarray(shift_tmp);
    }
    if (PyIs_NativeInt(nbproc_tmp))
    {
        nbproc = PyInt64_to_Int64(nbproc_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    if (!pyarray_check(w_x_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_x = pyarray_to_ndarray(w_x_tmp);
    }
    if (!pyarray_check(w_y_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_y = pyarray_to_ndarray(w_y_tmp);
    }
    if (!pyarray_check(w_z_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_z = pyarray_to_ndarray(w_z_tmp);
    }
    bind_c_cell_gradient_2d(nd_ndim(&w_c, 0), nd_data(&w_c), nd_ndim(&w_ghost, 0), nd_data(&w_ghost), nd_ndim(&w_halo, 0), nd_data(&w_halo), nd_ndim(&w_haloghost, 0), nd_data(&w_haloghost), nd_ndim(&centerc, 0), nd_ndim(&centerc, 1), nd_data(&centerc), nd_ndim(&cellnid, 0), nd_ndim(&cellnid, 1), nd_data(&cellnid), nd_ndim(&halonid, 0), nd_ndim(&halonid, 1), nd_data(&halonid), nd_ndim(&nodecid, 0), nd_ndim(&nodecid, 1), nd_data(&nodecid), nd_ndim(&periodicn, 0), nd_ndim(&periodicn, 1), nd_data(&periodicn), nd_ndim(&periodic, 0), nd_ndim(&periodic, 1), nd_data(&periodic), nd_ndim(&namen, 0), nd_data(&namen), nd_ndim(&centerg, 0), nd_ndim(&centerg, 1), nd_ndim(&centerg, 2), nd_data(&centerg), nd_ndim(&halocenterg, 0), nd_ndim(&halocenterg, 1), nd_ndim(&halocenterg, 2), nd_data(&halocenterg), nd_ndim(&vertexn, 0), nd_ndim(&vertexn, 1), nd_data(&vertexn), nd_ndim(&centerh, 0), nd_ndim(&centerh, 1), nd_data(&centerh), nd_ndim(&shift, 0), nd_ndim(&shift, 1), nd_data(&shift), nbproc, nd_ndim(&w_x, 0), nd_data(&w_x), nd_ndim(&w_y, 0), nd_data(&w_y), nd_ndim(&w_z, 0), nd_data(&w_z));
    result = Py_BuildValue("");
    free_pointer(w_c);
    free_pointer(w_ghost);
    free_pointer(w_halo);
    free_pointer(w_haloghost);
    free_pointer(centerc);
    free_pointer(cellnid);
    free_pointer(halonid);
    free_pointer(nodecid);
    free_pointer(periodicn);
    free_pointer(periodic);
    free_pointer(namen);
    free_pointer(centerg);
    free_pointer(halocenterg);
    free_pointer(vertexn);
    free_pointer(centerh);
    free_pointer(shift);
    free_pointer(w_x);
    free_pointer(w_y);
    free_pointer(w_z);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *cell_gradient_3d_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray w_c = {.shape = NULL};
    t_ndarray w_ghost = {.shape = NULL};
    t_ndarray w_halo = {.shape = NULL};
    t_ndarray w_haloghost = {.shape = NULL};
    t_ndarray centerc = {.shape = NULL};
    t_ndarray cellnid = {.shape = NULL};
    t_ndarray halonid = {.shape = NULL};
    t_ndarray nodecid = {.shape = NULL};
    t_ndarray periodicn = {.shape = NULL};
    t_ndarray periodic = {.shape = NULL};
    t_ndarray namen = {.shape = NULL};
    t_ndarray centerg = {.shape = NULL};
    t_ndarray halocenterg = {.shape = NULL};
    t_ndarray vertexn = {.shape = NULL};
    t_ndarray centerh = {.shape = NULL};
    t_ndarray shift = {.shape = NULL};
    int64_t nbproc;
    t_ndarray w_x = {.shape = NULL};
    t_ndarray w_y = {.shape = NULL};
    t_ndarray w_z = {.shape = NULL};
    PyArrayObject *w_c_tmp;
    PyArrayObject *w_ghost_tmp;
    PyArrayObject *w_halo_tmp;
    PyArrayObject *w_haloghost_tmp;
    PyArrayObject *centerc_tmp;
    PyArrayObject *cellnid_tmp;
    PyArrayObject *halonid_tmp;
    PyArrayObject *nodecid_tmp;
    PyArrayObject *periodicn_tmp;
    PyArrayObject *periodic_tmp;
    PyArrayObject *namen_tmp;
    PyArrayObject *centerg_tmp;
    PyArrayObject *halocenterg_tmp;
    PyArrayObject *vertexn_tmp;
    PyArrayObject *centerh_tmp;
    PyArrayObject *shift_tmp;
    PyObject *nbproc_tmp;
    PyArrayObject *w_x_tmp;
    PyArrayObject *w_y_tmp;
    PyArrayObject *w_z_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "w_c",
        "w_ghost",
        "w_halo",
        "w_haloghost",
        "centerc",
        "cellnid",
        "halonid",
        "nodecid",
        "periodicn",
        "periodic",
        "namen",
        "centerg",
        "halocenterg",
        "vertexn",
        "centerh",
        "shift",
        "nbproc",
        "w_x",
        "w_y",
        "w_z",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!OO!O!O!", kwlist, &PyArray_Type, &w_c_tmp, &PyArray_Type, &w_ghost_tmp, &PyArray_Type, &w_halo_tmp, &PyArray_Type, &w_haloghost_tmp, &PyArray_Type, &centerc_tmp, &PyArray_Type, &cellnid_tmp, &PyArray_Type, &halonid_tmp, &PyArray_Type, &nodecid_tmp, &PyArray_Type, &periodicn_tmp, &PyArray_Type, &periodic_tmp, &PyArray_Type, &namen_tmp, &PyArray_Type, &centerg_tmp, &PyArray_Type, &halocenterg_tmp, &PyArray_Type, &vertexn_tmp, &PyArray_Type, &centerh_tmp, &PyArray_Type, &shift_tmp, &nbproc_tmp, &PyArray_Type, &w_x_tmp, &PyArray_Type, &w_y_tmp, &PyArray_Type, &w_z_tmp))
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
    if (!pyarray_check(w_ghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_ghost = pyarray_to_ndarray(w_ghost_tmp);
    }
    if (!pyarray_check(w_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_halo = pyarray_to_ndarray(w_halo_tmp);
    }
    if (!pyarray_check(w_haloghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_haloghost = pyarray_to_ndarray(w_haloghost_tmp);
    }
    if (!pyarray_check(centerc_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerc = pyarray_to_ndarray(centerc_tmp);
    }
    if (!pyarray_check(cellnid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellnid = pyarray_to_ndarray(cellnid_tmp);
    }
    if (!pyarray_check(halonid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        halonid = pyarray_to_ndarray(halonid_tmp);
    }
    if (!pyarray_check(nodecid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nodecid = pyarray_to_ndarray(nodecid_tmp);
    }
    if (!pyarray_check(periodicn_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        periodicn = pyarray_to_ndarray(periodicn_tmp);
    }
    if (!pyarray_check(periodic_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        periodic = pyarray_to_ndarray(periodic_tmp);
    }
    if (!pyarray_check(namen_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        namen = pyarray_to_ndarray(namen_tmp);
    }
    if (!pyarray_check(centerg_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerg = pyarray_to_ndarray(centerg_tmp);
    }
    if (!pyarray_check(halocenterg_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        halocenterg = pyarray_to_ndarray(halocenterg_tmp);
    }
    if (!pyarray_check(vertexn_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        vertexn = pyarray_to_ndarray(vertexn_tmp);
    }
    if (!pyarray_check(centerh_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerh = pyarray_to_ndarray(centerh_tmp);
    }
    if (!pyarray_check(shift_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        shift = pyarray_to_ndarray(shift_tmp);
    }
    if (PyIs_NativeInt(nbproc_tmp))
    {
        nbproc = PyInt64_to_Int64(nbproc_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    if (!pyarray_check(w_x_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_x = pyarray_to_ndarray(w_x_tmp);
    }
    if (!pyarray_check(w_y_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_y = pyarray_to_ndarray(w_y_tmp);
    }
    if (!pyarray_check(w_z_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_z = pyarray_to_ndarray(w_z_tmp);
    }
    bind_c_cell_gradient_3d(nd_ndim(&w_c, 0), nd_data(&w_c), nd_ndim(&w_ghost, 0), nd_data(&w_ghost), nd_ndim(&w_halo, 0), nd_data(&w_halo), nd_ndim(&w_haloghost, 0), nd_data(&w_haloghost), nd_ndim(&centerc, 0), nd_ndim(&centerc, 1), nd_data(&centerc), nd_ndim(&cellnid, 0), nd_ndim(&cellnid, 1), nd_data(&cellnid), nd_ndim(&halonid, 0), nd_ndim(&halonid, 1), nd_data(&halonid), nd_ndim(&nodecid, 0), nd_ndim(&nodecid, 1), nd_data(&nodecid), nd_ndim(&periodicn, 0), nd_ndim(&periodicn, 1), nd_data(&periodicn), nd_ndim(&periodic, 0), nd_ndim(&periodic, 1), nd_data(&periodic), nd_ndim(&namen, 0), nd_data(&namen), nd_ndim(&centerg, 0), nd_ndim(&centerg, 1), nd_ndim(&centerg, 2), nd_data(&centerg), nd_ndim(&halocenterg, 0), nd_ndim(&halocenterg, 1), nd_ndim(&halocenterg, 2), nd_data(&halocenterg), nd_ndim(&vertexn, 0), nd_ndim(&vertexn, 1), nd_data(&vertexn), nd_ndim(&centerh, 0), nd_ndim(&centerh, 1), nd_data(&centerh), nd_ndim(&shift, 0), nd_ndim(&shift, 1), nd_data(&shift), nbproc, nd_ndim(&w_x, 0), nd_data(&w_x), nd_ndim(&w_y, 0), nd_data(&w_y), nd_ndim(&w_z, 0), nd_data(&w_z));
    result = Py_BuildValue("");
    free_pointer(w_c);
    free_pointer(w_ghost);
    free_pointer(w_halo);
    free_pointer(w_haloghost);
    free_pointer(centerc);
    free_pointer(cellnid);
    free_pointer(halonid);
    free_pointer(nodecid);
    free_pointer(periodicn);
    free_pointer(periodic);
    free_pointer(namen);
    free_pointer(centerg);
    free_pointer(halocenterg);
    free_pointer(vertexn);
    free_pointer(centerh);
    free_pointer(shift);
    free_pointer(w_x);
    free_pointer(w_y);
    free_pointer(w_z);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *face_gradient_2d_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray w_c = {.shape = NULL};
    t_ndarray w_ghost = {.shape = NULL};
    t_ndarray w_halo = {.shape = NULL};
    t_ndarray w_node = {.shape = NULL};
    t_ndarray cellidf = {.shape = NULL};
    t_ndarray nodeidf = {.shape = NULL};
    t_ndarray centergf = {.shape = NULL};
    t_ndarray namef = {.shape = NULL};
    t_ndarray halofid = {.shape = NULL};
    t_ndarray centerc = {.shape = NULL};
    t_ndarray centerh = {.shape = NULL};
    t_ndarray vertexn = {.shape = NULL};
    t_ndarray airDiamond = {.shape = NULL};
    t_ndarray normalf = {.shape = NULL};
    t_ndarray f_1 = {.shape = NULL};
    t_ndarray f_2 = {.shape = NULL};
    t_ndarray f_3 = {.shape = NULL};
    t_ndarray f_4 = {.shape = NULL};
    t_ndarray shift = {.shape = NULL};
    t_ndarray wx_face = {.shape = NULL};
    t_ndarray wy_face = {.shape = NULL};
    t_ndarray wz_face = {.shape = NULL};
    t_ndarray innerfaces = {.shape = NULL};
    t_ndarray halofaces = {.shape = NULL};
    t_ndarray dirichletfaces = {.shape = NULL};
    t_ndarray neumann = {.shape = NULL};
    t_ndarray periodicfaces = {.shape = NULL};
    PyArrayObject *w_c_tmp;
    PyArrayObject *w_ghost_tmp;
    PyArrayObject *w_halo_tmp;
    PyArrayObject *w_node_tmp;
    PyArrayObject *cellidf_tmp;
    PyArrayObject *nodeidf_tmp;
    PyArrayObject *centergf_tmp;
    PyArrayObject *namef_tmp;
    PyArrayObject *halofid_tmp;
    PyArrayObject *centerc_tmp;
    PyArrayObject *centerh_tmp;
    PyArrayObject *vertexn_tmp;
    PyArrayObject *airDiamond_tmp;
    PyArrayObject *normalf_tmp;
    PyArrayObject *f_1_tmp;
    PyArrayObject *f_2_tmp;
    PyArrayObject *f_3_tmp;
    PyArrayObject *f_4_tmp;
    PyArrayObject *shift_tmp;
    PyArrayObject *wx_face_tmp;
    PyArrayObject *wy_face_tmp;
    PyArrayObject *wz_face_tmp;
    PyArrayObject *innerfaces_tmp;
    PyArrayObject *halofaces_tmp;
    PyArrayObject *dirichletfaces_tmp;
    PyArrayObject *neumann_tmp;
    PyArrayObject *periodicfaces_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "w_c",
        "w_ghost",
        "w_halo",
        "w_node",
        "cellidf",
        "nodeidf",
        "centergf",
        "namef",
        "halofid",
        "centerc",
        "centerh",
        "vertexn",
        "airDiamond",
        "normalf",
        "f_1",
        "f_2",
        "f_3",
        "f_4",
        "shift",
        "wx_face",
        "wy_face",
        "wz_face",
        "innerfaces",
        "halofaces",
        "dirichletfaces",
        "neumann",
        "periodicfaces",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!", kwlist, &PyArray_Type, &w_c_tmp, &PyArray_Type, &w_ghost_tmp, &PyArray_Type, &w_halo_tmp, &PyArray_Type, &w_node_tmp, &PyArray_Type, &cellidf_tmp, &PyArray_Type, &nodeidf_tmp, &PyArray_Type, &centergf_tmp, &PyArray_Type, &namef_tmp, &PyArray_Type, &halofid_tmp, &PyArray_Type, &centerc_tmp, &PyArray_Type, &centerh_tmp, &PyArray_Type, &vertexn_tmp, &PyArray_Type, &airDiamond_tmp, &PyArray_Type, &normalf_tmp, &PyArray_Type, &f_1_tmp, &PyArray_Type, &f_2_tmp, &PyArray_Type, &f_3_tmp, &PyArray_Type, &f_4_tmp, &PyArray_Type, &shift_tmp, &PyArray_Type, &wx_face_tmp, &PyArray_Type, &wy_face_tmp, &PyArray_Type, &wz_face_tmp, &PyArray_Type, &innerfaces_tmp, &PyArray_Type, &halofaces_tmp, &PyArray_Type, &dirichletfaces_tmp, &PyArray_Type, &neumann_tmp, &PyArray_Type, &periodicfaces_tmp))
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
    if (!pyarray_check(w_ghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_ghost = pyarray_to_ndarray(w_ghost_tmp);
    }
    if (!pyarray_check(w_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_halo = pyarray_to_ndarray(w_halo_tmp);
    }
    if (!pyarray_check(w_node_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_node = pyarray_to_ndarray(w_node_tmp);
    }
    if (!pyarray_check(cellidf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellidf = pyarray_to_ndarray(cellidf_tmp);
    }
    if (!pyarray_check(nodeidf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nodeidf = pyarray_to_ndarray(nodeidf_tmp);
    }
    if (!pyarray_check(centergf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centergf = pyarray_to_ndarray(centergf_tmp);
    }
    if (!pyarray_check(namef_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        namef = pyarray_to_ndarray(namef_tmp);
    }
    if (!pyarray_check(halofid_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofid = pyarray_to_ndarray(halofid_tmp);
    }
    if (!pyarray_check(centerc_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerc = pyarray_to_ndarray(centerc_tmp);
    }
    if (!pyarray_check(centerh_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerh = pyarray_to_ndarray(centerh_tmp);
    }
    if (!pyarray_check(vertexn_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        vertexn = pyarray_to_ndarray(vertexn_tmp);
    }
    if (!pyarray_check(airDiamond_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        airDiamond = pyarray_to_ndarray(airDiamond_tmp);
    }
    if (!pyarray_check(normalf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        normalf = pyarray_to_ndarray(normalf_tmp);
    }
    if (!pyarray_check(f_1_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        f_1 = pyarray_to_ndarray(f_1_tmp);
    }
    if (!pyarray_check(f_2_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        f_2 = pyarray_to_ndarray(f_2_tmp);
    }
    if (!pyarray_check(f_3_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        f_3 = pyarray_to_ndarray(f_3_tmp);
    }
    if (!pyarray_check(f_4_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        f_4 = pyarray_to_ndarray(f_4_tmp);
    }
    if (!pyarray_check(shift_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        shift = pyarray_to_ndarray(shift_tmp);
    }
    if (!pyarray_check(wx_face_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        wx_face = pyarray_to_ndarray(wx_face_tmp);
    }
    if (!pyarray_check(wy_face_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        wy_face = pyarray_to_ndarray(wy_face_tmp);
    }
    if (!pyarray_check(wz_face_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        wz_face = pyarray_to_ndarray(wz_face_tmp);
    }
    if (!pyarray_check(innerfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        innerfaces = pyarray_to_ndarray(innerfaces_tmp);
    }
    if (!pyarray_check(halofaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofaces = pyarray_to_ndarray(halofaces_tmp);
    }
    if (!pyarray_check(dirichletfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        dirichletfaces = pyarray_to_ndarray(dirichletfaces_tmp);
    }
    if (!pyarray_check(neumann_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        neumann = pyarray_to_ndarray(neumann_tmp);
    }
    if (!pyarray_check(periodicfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        periodicfaces = pyarray_to_ndarray(periodicfaces_tmp);
    }
    bind_c_face_gradient_2d(nd_ndim(&w_c, 0), nd_data(&w_c), nd_ndim(&w_ghost, 0), nd_data(&w_ghost), nd_ndim(&w_halo, 0), nd_data(&w_halo), nd_ndim(&w_node, 0), nd_data(&w_node), nd_ndim(&cellidf, 0), nd_ndim(&cellidf, 1), nd_data(&cellidf), nd_ndim(&nodeidf, 0), nd_ndim(&nodeidf, 1), nd_data(&nodeidf), nd_ndim(&centergf, 0), nd_ndim(&centergf, 1), nd_data(&centergf), nd_ndim(&namef, 0), nd_data(&namef), nd_ndim(&halofid, 0), nd_data(&halofid), nd_ndim(&centerc, 0), nd_ndim(&centerc, 1), nd_data(&centerc), nd_ndim(&centerh, 0), nd_ndim(&centerh, 1), nd_data(&centerh), nd_ndim(&vertexn, 0), nd_ndim(&vertexn, 1), nd_data(&vertexn), nd_ndim(&airDiamond, 0), nd_data(&airDiamond), nd_ndim(&normalf, 0), nd_ndim(&normalf, 1), nd_data(&normalf), nd_ndim(&f_1, 0), nd_ndim(&f_1, 1), nd_data(&f_1), nd_ndim(&f_2, 0), nd_ndim(&f_2, 1), nd_data(&f_2), nd_ndim(&f_3, 0), nd_ndim(&f_3, 1), nd_data(&f_3), nd_ndim(&f_4, 0), nd_ndim(&f_4, 1), nd_data(&f_4), nd_ndim(&shift, 0), nd_ndim(&shift, 1), nd_data(&shift), nd_ndim(&wx_face, 0), nd_data(&wx_face), nd_ndim(&wy_face, 0), nd_data(&wy_face), nd_ndim(&wz_face, 0), nd_data(&wz_face), nd_ndim(&innerfaces, 0), nd_data(&innerfaces), nd_ndim(&halofaces, 0), nd_data(&halofaces), nd_ndim(&dirichletfaces, 0), nd_data(&dirichletfaces), nd_ndim(&neumann, 0), nd_data(&neumann), nd_ndim(&periodicfaces, 0), nd_data(&periodicfaces));
    result = Py_BuildValue("");
    free_pointer(w_c);
    free_pointer(w_ghost);
    free_pointer(w_halo);
    free_pointer(w_node);
    free_pointer(cellidf);
    free_pointer(nodeidf);
    free_pointer(centergf);
    free_pointer(namef);
    free_pointer(halofid);
    free_pointer(centerc);
    free_pointer(centerh);
    free_pointer(vertexn);
    free_pointer(airDiamond);
    free_pointer(normalf);
    free_pointer(f_1);
    free_pointer(f_2);
    free_pointer(f_3);
    free_pointer(f_4);
    free_pointer(shift);
    free_pointer(wx_face);
    free_pointer(wy_face);
    free_pointer(wz_face);
    free_pointer(innerfaces);
    free_pointer(halofaces);
    free_pointer(dirichletfaces);
    free_pointer(neumann);
    free_pointer(periodicfaces);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *face_gradient_2d_uv_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray h_c = {.shape = NULL};
    t_ndarray h_ghost = {.shape = NULL};
    t_ndarray h_halo = {.shape = NULL};
    t_ndarray h_node = {.shape = NULL};
    t_ndarray w_c = {.shape = NULL};
    t_ndarray w_ghost = {.shape = NULL};
    t_ndarray w_halo = {.shape = NULL};
    t_ndarray w_node = {.shape = NULL};
    t_ndarray cellidf = {.shape = NULL};
    t_ndarray nodeidf = {.shape = NULL};
    t_ndarray centergf = {.shape = NULL};
    t_ndarray namef = {.shape = NULL};
    t_ndarray halofid = {.shape = NULL};
    t_ndarray centerc = {.shape = NULL};
    t_ndarray centerh = {.shape = NULL};
    t_ndarray vertexn = {.shape = NULL};
    t_ndarray airDiamond = {.shape = NULL};
    t_ndarray normalf = {.shape = NULL};
    t_ndarray f_1 = {.shape = NULL};
    t_ndarray f_2 = {.shape = NULL};
    t_ndarray f_3 = {.shape = NULL};
    t_ndarray f_4 = {.shape = NULL};
    t_ndarray shift = {.shape = NULL};
    t_ndarray wx_face = {.shape = NULL};
    t_ndarray wy_face = {.shape = NULL};
    t_ndarray wz_face = {.shape = NULL};
    t_ndarray innerfaces = {.shape = NULL};
    t_ndarray halofaces = {.shape = NULL};
    t_ndarray dirichletfaces = {.shape = NULL};
    t_ndarray neumann = {.shape = NULL};
    t_ndarray periodicfaces = {.shape = NULL};
    PyArrayObject *h_c_tmp;
    PyArrayObject *h_ghost_tmp;
    PyArrayObject *h_halo_tmp;
    PyArrayObject *h_node_tmp;
    PyArrayObject *w_c_tmp;
    PyArrayObject *w_ghost_tmp;
    PyArrayObject *w_halo_tmp;
    PyArrayObject *w_node_tmp;
    PyArrayObject *cellidf_tmp;
    PyArrayObject *nodeidf_tmp;
    PyArrayObject *centergf_tmp;
    PyArrayObject *namef_tmp;
    PyArrayObject *halofid_tmp;
    PyArrayObject *centerc_tmp;
    PyArrayObject *centerh_tmp;
    PyArrayObject *vertexn_tmp;
    PyArrayObject *airDiamond_tmp;
    PyArrayObject *normalf_tmp;
    PyArrayObject *f_1_tmp;
    PyArrayObject *f_2_tmp;
    PyArrayObject *f_3_tmp;
    PyArrayObject *f_4_tmp;
    PyArrayObject *shift_tmp;
    PyArrayObject *wx_face_tmp;
    PyArrayObject *wy_face_tmp;
    PyArrayObject *wz_face_tmp;
    PyArrayObject *innerfaces_tmp;
    PyArrayObject *halofaces_tmp;
    PyArrayObject *dirichletfaces_tmp;
    PyArrayObject *neumann_tmp;
    PyArrayObject *periodicfaces_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "h_c",
        "h_ghost",
        "h_halo",
        "h_node",
        "w_c",
        "w_ghost",
        "w_halo",
        "w_node",
        "cellidf",
        "nodeidf",
        "centergf",
        "namef",
        "halofid",
        "centerc",
        "centerh",
        "vertexn",
        "airDiamond",
        "normalf",
        "f_1",
        "f_2",
        "f_3",
        "f_4",
        "shift",
        "wx_face",
        "wy_face",
        "wz_face",
        "innerfaces",
        "halofaces",
        "dirichletfaces",
        "neumann",
        "periodicfaces",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!", kwlist, &PyArray_Type, &h_c_tmp, &PyArray_Type, &h_ghost_tmp, &PyArray_Type, &h_halo_tmp, &PyArray_Type, &h_node_tmp, &PyArray_Type, &w_c_tmp, &PyArray_Type, &w_ghost_tmp, &PyArray_Type, &w_halo_tmp, &PyArray_Type, &w_node_tmp, &PyArray_Type, &cellidf_tmp, &PyArray_Type, &nodeidf_tmp, &PyArray_Type, &centergf_tmp, &PyArray_Type, &namef_tmp, &PyArray_Type, &halofid_tmp, &PyArray_Type, &centerc_tmp, &PyArray_Type, &centerh_tmp, &PyArray_Type, &vertexn_tmp, &PyArray_Type, &airDiamond_tmp, &PyArray_Type, &normalf_tmp, &PyArray_Type, &f_1_tmp, &PyArray_Type, &f_2_tmp, &PyArray_Type, &f_3_tmp, &PyArray_Type, &f_4_tmp, &PyArray_Type, &shift_tmp, &PyArray_Type, &wx_face_tmp, &PyArray_Type, &wy_face_tmp, &PyArray_Type, &wz_face_tmp, &PyArray_Type, &innerfaces_tmp, &PyArray_Type, &halofaces_tmp, &PyArray_Type, &dirichletfaces_tmp, &PyArray_Type, &neumann_tmp, &PyArray_Type, &periodicfaces_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(h_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        h_c = pyarray_to_ndarray(h_c_tmp);
    }
    if (!pyarray_check(h_ghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        h_ghost = pyarray_to_ndarray(h_ghost_tmp);
    }
    if (!pyarray_check(h_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        h_halo = pyarray_to_ndarray(h_halo_tmp);
    }
    if (!pyarray_check(h_node_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        h_node = pyarray_to_ndarray(h_node_tmp);
    }
    if (!pyarray_check(w_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_c = pyarray_to_ndarray(w_c_tmp);
    }
    if (!pyarray_check(w_ghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_ghost = pyarray_to_ndarray(w_ghost_tmp);
    }
    if (!pyarray_check(w_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_halo = pyarray_to_ndarray(w_halo_tmp);
    }
    if (!pyarray_check(w_node_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_node = pyarray_to_ndarray(w_node_tmp);
    }
    if (!pyarray_check(cellidf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellidf = pyarray_to_ndarray(cellidf_tmp);
    }
    if (!pyarray_check(nodeidf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nodeidf = pyarray_to_ndarray(nodeidf_tmp);
    }
    if (!pyarray_check(centergf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centergf = pyarray_to_ndarray(centergf_tmp);
    }
    if (!pyarray_check(namef_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        namef = pyarray_to_ndarray(namef_tmp);
    }
    if (!pyarray_check(halofid_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofid = pyarray_to_ndarray(halofid_tmp);
    }
    if (!pyarray_check(centerc_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerc = pyarray_to_ndarray(centerc_tmp);
    }
    if (!pyarray_check(centerh_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerh = pyarray_to_ndarray(centerh_tmp);
    }
    if (!pyarray_check(vertexn_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        vertexn = pyarray_to_ndarray(vertexn_tmp);
    }
    if (!pyarray_check(airDiamond_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        airDiamond = pyarray_to_ndarray(airDiamond_tmp);
    }
    if (!pyarray_check(normalf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        normalf = pyarray_to_ndarray(normalf_tmp);
    }
    if (!pyarray_check(f_1_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        f_1 = pyarray_to_ndarray(f_1_tmp);
    }
    if (!pyarray_check(f_2_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        f_2 = pyarray_to_ndarray(f_2_tmp);
    }
    if (!pyarray_check(f_3_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        f_3 = pyarray_to_ndarray(f_3_tmp);
    }
    if (!pyarray_check(f_4_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        f_4 = pyarray_to_ndarray(f_4_tmp);
    }
    if (!pyarray_check(shift_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        shift = pyarray_to_ndarray(shift_tmp);
    }
    if (!pyarray_check(wx_face_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        wx_face = pyarray_to_ndarray(wx_face_tmp);
    }
    if (!pyarray_check(wy_face_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        wy_face = pyarray_to_ndarray(wy_face_tmp);
    }
    if (!pyarray_check(wz_face_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        wz_face = pyarray_to_ndarray(wz_face_tmp);
    }
    if (!pyarray_check(innerfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        innerfaces = pyarray_to_ndarray(innerfaces_tmp);
    }
    if (!pyarray_check(halofaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofaces = pyarray_to_ndarray(halofaces_tmp);
    }
    if (!pyarray_check(dirichletfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        dirichletfaces = pyarray_to_ndarray(dirichletfaces_tmp);
    }
    if (!pyarray_check(neumann_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        neumann = pyarray_to_ndarray(neumann_tmp);
    }
    if (!pyarray_check(periodicfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        periodicfaces = pyarray_to_ndarray(periodicfaces_tmp);
    }
    bind_c_face_gradient_2d_uv(nd_ndim(&h_c, 0), nd_data(&h_c), nd_ndim(&h_ghost, 0), nd_data(&h_ghost), nd_ndim(&h_halo, 0), nd_data(&h_halo), nd_ndim(&h_node, 0), nd_data(&h_node), nd_ndim(&w_c, 0), nd_data(&w_c), nd_ndim(&w_ghost, 0), nd_data(&w_ghost), nd_ndim(&w_halo, 0), nd_data(&w_halo), nd_ndim(&w_node, 0), nd_data(&w_node), nd_ndim(&cellidf, 0), nd_ndim(&cellidf, 1), nd_data(&cellidf), nd_ndim(&nodeidf, 0), nd_ndim(&nodeidf, 1), nd_data(&nodeidf), nd_ndim(&centergf, 0), nd_ndim(&centergf, 1), nd_data(&centergf), nd_ndim(&namef, 0), nd_data(&namef), nd_ndim(&halofid, 0), nd_data(&halofid), nd_ndim(&centerc, 0), nd_ndim(&centerc, 1), nd_data(&centerc), nd_ndim(&centerh, 0), nd_ndim(&centerh, 1), nd_data(&centerh), nd_ndim(&vertexn, 0), nd_ndim(&vertexn, 1), nd_data(&vertexn), nd_ndim(&airDiamond, 0), nd_data(&airDiamond), nd_ndim(&normalf, 0), nd_ndim(&normalf, 1), nd_data(&normalf), nd_ndim(&f_1, 0), nd_ndim(&f_1, 1), nd_data(&f_1), nd_ndim(&f_2, 0), nd_ndim(&f_2, 1), nd_data(&f_2), nd_ndim(&f_3, 0), nd_ndim(&f_3, 1), nd_data(&f_3), nd_ndim(&f_4, 0), nd_ndim(&f_4, 1), nd_data(&f_4), nd_ndim(&shift, 0), nd_ndim(&shift, 1), nd_data(&shift), nd_ndim(&wx_face, 0), nd_data(&wx_face), nd_ndim(&wy_face, 0), nd_data(&wy_face), nd_ndim(&wz_face, 0), nd_data(&wz_face), nd_ndim(&innerfaces, 0), nd_data(&innerfaces), nd_ndim(&halofaces, 0), nd_data(&halofaces), nd_ndim(&dirichletfaces, 0), nd_data(&dirichletfaces), nd_ndim(&neumann, 0), nd_data(&neumann), nd_ndim(&periodicfaces, 0), nd_data(&periodicfaces));
    result = Py_BuildValue("");
    free_pointer(h_c);
    free_pointer(h_ghost);
    free_pointer(h_halo);
    free_pointer(h_node);
    free_pointer(w_c);
    free_pointer(w_ghost);
    free_pointer(w_halo);
    free_pointer(w_node);
    free_pointer(cellidf);
    free_pointer(nodeidf);
    free_pointer(centergf);
    free_pointer(namef);
    free_pointer(halofid);
    free_pointer(centerc);
    free_pointer(centerh);
    free_pointer(vertexn);
    free_pointer(airDiamond);
    free_pointer(normalf);
    free_pointer(f_1);
    free_pointer(f_2);
    free_pointer(f_3);
    free_pointer(f_4);
    free_pointer(shift);
    free_pointer(wx_face);
    free_pointer(wy_face);
    free_pointer(wz_face);
    free_pointer(innerfaces);
    free_pointer(halofaces);
    free_pointer(dirichletfaces);
    free_pointer(neumann);
    free_pointer(periodicfaces);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *face_gradient_3d_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray w_c = {.shape = NULL};
    t_ndarray w_ghost = {.shape = NULL};
    t_ndarray w_halo = {.shape = NULL};
    t_ndarray w_node = {.shape = NULL};
    t_ndarray cellidf = {.shape = NULL};
    t_ndarray nodeidf = {.shape = NULL};
    t_ndarray centergf = {.shape = NULL};
    t_ndarray namef = {.shape = NULL};
    t_ndarray halofid = {.shape = NULL};
    t_ndarray centerc = {.shape = NULL};
    t_ndarray centerh = {.shape = NULL};
    t_ndarray vertexn = {.shape = NULL};
    t_ndarray airDiamond = {.shape = NULL};
    t_ndarray normalf = {.shape = NULL};
    t_ndarray f_1 = {.shape = NULL};
    t_ndarray f_2 = {.shape = NULL};
    t_ndarray f_3 = {.shape = NULL};
    t_ndarray f_4 = {.shape = NULL};
    t_ndarray shift = {.shape = NULL};
    t_ndarray wx_face = {.shape = NULL};
    t_ndarray wy_face = {.shape = NULL};
    t_ndarray wz_face = {.shape = NULL};
    t_ndarray innerfaces = {.shape = NULL};
    t_ndarray halofaces = {.shape = NULL};
    t_ndarray dirichletfaces = {.shape = NULL};
    t_ndarray neumann = {.shape = NULL};
    t_ndarray periodicfaces = {.shape = NULL};
    PyArrayObject *w_c_tmp;
    PyArrayObject *w_ghost_tmp;
    PyArrayObject *w_halo_tmp;
    PyArrayObject *w_node_tmp;
    PyArrayObject *cellidf_tmp;
    PyArrayObject *nodeidf_tmp;
    PyArrayObject *centergf_tmp;
    PyArrayObject *namef_tmp;
    PyArrayObject *halofid_tmp;
    PyArrayObject *centerc_tmp;
    PyArrayObject *centerh_tmp;
    PyArrayObject *vertexn_tmp;
    PyArrayObject *airDiamond_tmp;
    PyArrayObject *normalf_tmp;
    PyArrayObject *f_1_tmp;
    PyArrayObject *f_2_tmp;
    PyArrayObject *f_3_tmp;
    PyArrayObject *f_4_tmp;
    PyArrayObject *shift_tmp;
    PyArrayObject *wx_face_tmp;
    PyArrayObject *wy_face_tmp;
    PyArrayObject *wz_face_tmp;
    PyArrayObject *innerfaces_tmp;
    PyArrayObject *halofaces_tmp;
    PyArrayObject *dirichletfaces_tmp;
    PyArrayObject *neumann_tmp;
    PyArrayObject *periodicfaces_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "w_c",
        "w_ghost",
        "w_halo",
        "w_node",
        "cellidf",
        "nodeidf",
        "centergf",
        "namef",
        "halofid",
        "centerc",
        "centerh",
        "vertexn",
        "airDiamond",
        "normalf",
        "f_1",
        "f_2",
        "f_3",
        "f_4",
        "shift",
        "wx_face",
        "wy_face",
        "wz_face",
        "innerfaces",
        "halofaces",
        "dirichletfaces",
        "neumann",
        "periodicfaces",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!", kwlist, &PyArray_Type, &w_c_tmp, &PyArray_Type, &w_ghost_tmp, &PyArray_Type, &w_halo_tmp, &PyArray_Type, &w_node_tmp, &PyArray_Type, &cellidf_tmp, &PyArray_Type, &nodeidf_tmp, &PyArray_Type, &centergf_tmp, &PyArray_Type, &namef_tmp, &PyArray_Type, &halofid_tmp, &PyArray_Type, &centerc_tmp, &PyArray_Type, &centerh_tmp, &PyArray_Type, &vertexn_tmp, &PyArray_Type, &airDiamond_tmp, &PyArray_Type, &normalf_tmp, &PyArray_Type, &f_1_tmp, &PyArray_Type, &f_2_tmp, &PyArray_Type, &f_3_tmp, &PyArray_Type, &f_4_tmp, &PyArray_Type, &shift_tmp, &PyArray_Type, &wx_face_tmp, &PyArray_Type, &wy_face_tmp, &PyArray_Type, &wz_face_tmp, &PyArray_Type, &innerfaces_tmp, &PyArray_Type, &halofaces_tmp, &PyArray_Type, &dirichletfaces_tmp, &PyArray_Type, &neumann_tmp, &PyArray_Type, &periodicfaces_tmp))
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
    if (!pyarray_check(w_ghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_ghost = pyarray_to_ndarray(w_ghost_tmp);
    }
    if (!pyarray_check(w_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_halo = pyarray_to_ndarray(w_halo_tmp);
    }
    if (!pyarray_check(w_node_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_node = pyarray_to_ndarray(w_node_tmp);
    }
    if (!pyarray_check(cellidf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellidf = pyarray_to_ndarray(cellidf_tmp);
    }
    if (!pyarray_check(nodeidf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nodeidf = pyarray_to_ndarray(nodeidf_tmp);
    }
    if (!pyarray_check(centergf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centergf = pyarray_to_ndarray(centergf_tmp);
    }
    if (!pyarray_check(namef_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        namef = pyarray_to_ndarray(namef_tmp);
    }
    if (!pyarray_check(halofid_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofid = pyarray_to_ndarray(halofid_tmp);
    }
    if (!pyarray_check(centerc_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerc = pyarray_to_ndarray(centerc_tmp);
    }
    if (!pyarray_check(centerh_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerh = pyarray_to_ndarray(centerh_tmp);
    }
    if (!pyarray_check(vertexn_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        vertexn = pyarray_to_ndarray(vertexn_tmp);
    }
    if (!pyarray_check(airDiamond_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        airDiamond = pyarray_to_ndarray(airDiamond_tmp);
    }
    if (!pyarray_check(normalf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        normalf = pyarray_to_ndarray(normalf_tmp);
    }
    if (!pyarray_check(f_1_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        f_1 = pyarray_to_ndarray(f_1_tmp);
    }
    if (!pyarray_check(f_2_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        f_2 = pyarray_to_ndarray(f_2_tmp);
    }
    if (!pyarray_check(f_3_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        f_3 = pyarray_to_ndarray(f_3_tmp);
    }
    if (!pyarray_check(f_4_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        f_4 = pyarray_to_ndarray(f_4_tmp);
    }
    if (!pyarray_check(shift_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        shift = pyarray_to_ndarray(shift_tmp);
    }
    if (!pyarray_check(wx_face_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        wx_face = pyarray_to_ndarray(wx_face_tmp);
    }
    if (!pyarray_check(wy_face_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        wy_face = pyarray_to_ndarray(wy_face_tmp);
    }
    if (!pyarray_check(wz_face_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        wz_face = pyarray_to_ndarray(wz_face_tmp);
    }
    if (!pyarray_check(innerfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        innerfaces = pyarray_to_ndarray(innerfaces_tmp);
    }
    if (!pyarray_check(halofaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofaces = pyarray_to_ndarray(halofaces_tmp);
    }
    if (!pyarray_check(dirichletfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        dirichletfaces = pyarray_to_ndarray(dirichletfaces_tmp);
    }
    if (!pyarray_check(neumann_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        neumann = pyarray_to_ndarray(neumann_tmp);
    }
    if (!pyarray_check(periodicfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        periodicfaces = pyarray_to_ndarray(periodicfaces_tmp);
    }
    bind_c_face_gradient_3d(nd_ndim(&w_c, 0), nd_data(&w_c), nd_ndim(&w_ghost, 0), nd_data(&w_ghost), nd_ndim(&w_halo, 0), nd_data(&w_halo), nd_ndim(&w_node, 0), nd_data(&w_node), nd_ndim(&cellidf, 0), nd_ndim(&cellidf, 1), nd_data(&cellidf), nd_ndim(&nodeidf, 0), nd_ndim(&nodeidf, 1), nd_data(&nodeidf), nd_ndim(&centergf, 0), nd_ndim(&centergf, 1), nd_data(&centergf), nd_ndim(&namef, 0), nd_data(&namef), nd_ndim(&halofid, 0), nd_data(&halofid), nd_ndim(&centerc, 0), nd_ndim(&centerc, 1), nd_data(&centerc), nd_ndim(&centerh, 0), nd_ndim(&centerh, 1), nd_data(&centerh), nd_ndim(&vertexn, 0), nd_ndim(&vertexn, 1), nd_data(&vertexn), nd_ndim(&airDiamond, 0), nd_data(&airDiamond), nd_ndim(&normalf, 0), nd_ndim(&normalf, 1), nd_data(&normalf), nd_ndim(&f_1, 0), nd_ndim(&f_1, 1), nd_data(&f_1), nd_ndim(&f_2, 0), nd_ndim(&f_2, 1), nd_data(&f_2), nd_ndim(&f_3, 0), nd_ndim(&f_3, 1), nd_data(&f_3), nd_ndim(&f_4, 0), nd_ndim(&f_4, 1), nd_data(&f_4), nd_ndim(&shift, 0), nd_ndim(&shift, 1), nd_data(&shift), nd_ndim(&wx_face, 0), nd_data(&wx_face), nd_ndim(&wy_face, 0), nd_data(&wy_face), nd_ndim(&wz_face, 0), nd_data(&wz_face), nd_ndim(&innerfaces, 0), nd_data(&innerfaces), nd_ndim(&halofaces, 0), nd_data(&halofaces), nd_ndim(&dirichletfaces, 0), nd_data(&dirichletfaces), nd_ndim(&neumann, 0), nd_data(&neumann), nd_ndim(&periodicfaces, 0), nd_data(&periodicfaces));
    result = Py_BuildValue("");
    free_pointer(w_c);
    free_pointer(w_ghost);
    free_pointer(w_halo);
    free_pointer(w_node);
    free_pointer(cellidf);
    free_pointer(nodeidf);
    free_pointer(centergf);
    free_pointer(namef);
    free_pointer(halofid);
    free_pointer(centerc);
    free_pointer(centerh);
    free_pointer(vertexn);
    free_pointer(airDiamond);
    free_pointer(normalf);
    free_pointer(f_1);
    free_pointer(f_2);
    free_pointer(f_3);
    free_pointer(f_4);
    free_pointer(shift);
    free_pointer(wx_face);
    free_pointer(wy_face);
    free_pointer(wz_face);
    free_pointer(innerfaces);
    free_pointer(halofaces);
    free_pointer(dirichletfaces);
    free_pointer(neumann);
    free_pointer(periodicfaces);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *centertovertex_2d_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray w_c = {.shape = NULL};
    t_ndarray w_ghost = {.shape = NULL};
    t_ndarray w_halo = {.shape = NULL};
    t_ndarray w_haloghost = {.shape = NULL};
    t_ndarray centerc = {.shape = NULL};
    t_ndarray centerh = {.shape = NULL};
    t_ndarray cellidn = {.shape = NULL};
    t_ndarray periodicn = {.shape = NULL};
    t_ndarray haloidn = {.shape = NULL};
    t_ndarray vertexn = {.shape = NULL};
    t_ndarray namen = {.shape = NULL};
    t_ndarray centergn = {.shape = NULL};
    t_ndarray halocentergn = {.shape = NULL};
    t_ndarray R_x = {.shape = NULL};
    t_ndarray R_y = {.shape = NULL};
    t_ndarray lambda_x = {.shape = NULL};
    t_ndarray lambda_y = {.shape = NULL};
    t_ndarray number = {.shape = NULL};
    t_ndarray shift = {.shape = NULL};
    int64_t nbproc;
    t_ndarray w_n = {.shape = NULL};
    PyArrayObject *w_c_tmp;
    PyArrayObject *w_ghost_tmp;
    PyArrayObject *w_halo_tmp;
    PyArrayObject *w_haloghost_tmp;
    PyArrayObject *centerc_tmp;
    PyArrayObject *centerh_tmp;
    PyArrayObject *cellidn_tmp;
    PyArrayObject *periodicn_tmp;
    PyArrayObject *haloidn_tmp;
    PyArrayObject *vertexn_tmp;
    PyArrayObject *namen_tmp;
    PyArrayObject *centergn_tmp;
    PyArrayObject *halocentergn_tmp;
    PyArrayObject *R_x_tmp;
    PyArrayObject *R_y_tmp;
    PyArrayObject *lambda_x_tmp;
    PyArrayObject *lambda_y_tmp;
    PyArrayObject *number_tmp;
    PyArrayObject *shift_tmp;
    PyObject *nbproc_tmp;
    PyArrayObject *w_n_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "w_c",
        "w_ghost",
        "w_halo",
        "w_haloghost",
        "centerc",
        "centerh",
        "cellidn",
        "periodicn",
        "haloidn",
        "vertexn",
        "namen",
        "centergn",
        "halocentergn",
        "R_x",
        "R_y",
        "lambda_x",
        "lambda_y",
        "number",
        "shift",
        "nbproc",
        "w_n",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!OO!", kwlist, &PyArray_Type, &w_c_tmp, &PyArray_Type, &w_ghost_tmp, &PyArray_Type, &w_halo_tmp, &PyArray_Type, &w_haloghost_tmp, &PyArray_Type, &centerc_tmp, &PyArray_Type, &centerh_tmp, &PyArray_Type, &cellidn_tmp, &PyArray_Type, &periodicn_tmp, &PyArray_Type, &haloidn_tmp, &PyArray_Type, &vertexn_tmp, &PyArray_Type, &namen_tmp, &PyArray_Type, &centergn_tmp, &PyArray_Type, &halocentergn_tmp, &PyArray_Type, &R_x_tmp, &PyArray_Type, &R_y_tmp, &PyArray_Type, &lambda_x_tmp, &PyArray_Type, &lambda_y_tmp, &PyArray_Type, &number_tmp, &PyArray_Type, &shift_tmp, &nbproc_tmp, &PyArray_Type, &w_n_tmp))
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
    if (!pyarray_check(w_ghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_ghost = pyarray_to_ndarray(w_ghost_tmp);
    }
    if (!pyarray_check(w_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_halo = pyarray_to_ndarray(w_halo_tmp);
    }
    if (!pyarray_check(w_haloghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_haloghost = pyarray_to_ndarray(w_haloghost_tmp);
    }
    if (!pyarray_check(centerc_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerc = pyarray_to_ndarray(centerc_tmp);
    }
    if (!pyarray_check(centerh_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerh = pyarray_to_ndarray(centerh_tmp);
    }
    if (!pyarray_check(cellidn_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellidn = pyarray_to_ndarray(cellidn_tmp);
    }
    if (!pyarray_check(periodicn_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        periodicn = pyarray_to_ndarray(periodicn_tmp);
    }
    if (!pyarray_check(haloidn_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        haloidn = pyarray_to_ndarray(haloidn_tmp);
    }
    if (!pyarray_check(vertexn_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        vertexn = pyarray_to_ndarray(vertexn_tmp);
    }
    if (!pyarray_check(namen_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        namen = pyarray_to_ndarray(namen_tmp);
    }
    if (!pyarray_check(centergn_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centergn = pyarray_to_ndarray(centergn_tmp);
    }
    if (!pyarray_check(halocentergn_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        halocentergn = pyarray_to_ndarray(halocentergn_tmp);
    }
    if (!pyarray_check(R_x_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        R_x = pyarray_to_ndarray(R_x_tmp);
    }
    if (!pyarray_check(R_y_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        R_y = pyarray_to_ndarray(R_y_tmp);
    }
    if (!pyarray_check(lambda_x_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        lambda_x = pyarray_to_ndarray(lambda_x_tmp);
    }
    if (!pyarray_check(lambda_y_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        lambda_y = pyarray_to_ndarray(lambda_y_tmp);
    }
    if (!pyarray_check(number_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        number = pyarray_to_ndarray(number_tmp);
    }
    if (!pyarray_check(shift_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        shift = pyarray_to_ndarray(shift_tmp);
    }
    if (PyIs_NativeInt(nbproc_tmp))
    {
        nbproc = PyInt64_to_Int64(nbproc_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    if (!pyarray_check(w_n_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_n = pyarray_to_ndarray(w_n_tmp);
    }
    bind_c_centertovertex_2d(nd_ndim(&w_c, 0), nd_data(&w_c), nd_ndim(&w_ghost, 0), nd_data(&w_ghost), nd_ndim(&w_halo, 0), nd_data(&w_halo), nd_ndim(&w_haloghost, 0), nd_data(&w_haloghost), nd_ndim(&centerc, 0), nd_ndim(&centerc, 1), nd_data(&centerc), nd_ndim(&centerh, 0), nd_ndim(&centerh, 1), nd_data(&centerh), nd_ndim(&cellidn, 0), nd_ndim(&cellidn, 1), nd_data(&cellidn), nd_ndim(&periodicn, 0), nd_ndim(&periodicn, 1), nd_data(&periodicn), nd_ndim(&haloidn, 0), nd_ndim(&haloidn, 1), nd_data(&haloidn), nd_ndim(&vertexn, 0), nd_ndim(&vertexn, 1), nd_data(&vertexn), nd_ndim(&namen, 0), nd_data(&namen), nd_ndim(&centergn, 0), nd_ndim(&centergn, 1), nd_ndim(&centergn, 2), nd_data(&centergn), nd_ndim(&halocentergn, 0), nd_ndim(&halocentergn, 1), nd_ndim(&halocentergn, 2), nd_data(&halocentergn), nd_ndim(&R_x, 0), nd_data(&R_x), nd_ndim(&R_y, 0), nd_data(&R_y), nd_ndim(&lambda_x, 0), nd_data(&lambda_x), nd_ndim(&lambda_y, 0), nd_data(&lambda_y), nd_ndim(&number, 0), nd_data(&number), nd_ndim(&shift, 0), nd_ndim(&shift, 1), nd_data(&shift), nbproc, nd_ndim(&w_n, 0), nd_data(&w_n));
    result = Py_BuildValue("");
    free_pointer(w_c);
    free_pointer(w_ghost);
    free_pointer(w_halo);
    free_pointer(w_haloghost);
    free_pointer(centerc);
    free_pointer(centerh);
    free_pointer(cellidn);
    free_pointer(periodicn);
    free_pointer(haloidn);
    free_pointer(vertexn);
    free_pointer(namen);
    free_pointer(centergn);
    free_pointer(halocentergn);
    free_pointer(R_x);
    free_pointer(R_y);
    free_pointer(lambda_x);
    free_pointer(lambda_y);
    free_pointer(number);
    free_pointer(shift);
    free_pointer(w_n);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *centertovertex_3d_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray w_c = {.shape = NULL};
    t_ndarray w_ghost = {.shape = NULL};
    t_ndarray w_halo = {.shape = NULL};
    t_ndarray w_haloghost = {.shape = NULL};
    t_ndarray centerc = {.shape = NULL};
    t_ndarray centerh = {.shape = NULL};
    t_ndarray cellidn = {.shape = NULL};
    t_ndarray periodicn = {.shape = NULL};
    t_ndarray haloidn = {.shape = NULL};
    t_ndarray vertexn = {.shape = NULL};
    t_ndarray namen = {.shape = NULL};
    t_ndarray centergn = {.shape = NULL};
    t_ndarray halocentergn = {.shape = NULL};
    t_ndarray R_x = {.shape = NULL};
    t_ndarray R_y = {.shape = NULL};
    t_ndarray R_z = {.shape = NULL};
    t_ndarray lambda_x = {.shape = NULL};
    t_ndarray lambda_y = {.shape = NULL};
    t_ndarray lambda_z = {.shape = NULL};
    t_ndarray number = {.shape = NULL};
    t_ndarray shift = {.shape = NULL};
    int64_t nbproc;
    t_ndarray w_n = {.shape = NULL};
    PyArrayObject *w_c_tmp;
    PyArrayObject *w_ghost_tmp;
    PyArrayObject *w_halo_tmp;
    PyArrayObject *w_haloghost_tmp;
    PyArrayObject *centerc_tmp;
    PyArrayObject *centerh_tmp;
    PyArrayObject *cellidn_tmp;
    PyArrayObject *periodicn_tmp;
    PyArrayObject *haloidn_tmp;
    PyArrayObject *vertexn_tmp;
    PyArrayObject *namen_tmp;
    PyArrayObject *centergn_tmp;
    PyArrayObject *halocentergn_tmp;
    PyArrayObject *R_x_tmp;
    PyArrayObject *R_y_tmp;
    PyArrayObject *R_z_tmp;
    PyArrayObject *lambda_x_tmp;
    PyArrayObject *lambda_y_tmp;
    PyArrayObject *lambda_z_tmp;
    PyArrayObject *number_tmp;
    PyArrayObject *shift_tmp;
    PyObject *nbproc_tmp;
    PyArrayObject *w_n_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "w_c",
        "w_ghost",
        "w_halo",
        "w_haloghost",
        "centerc",
        "centerh",
        "cellidn",
        "periodicn",
        "haloidn",
        "vertexn",
        "namen",
        "centergn",
        "halocentergn",
        "R_x",
        "R_y",
        "R_z",
        "lambda_x",
        "lambda_y",
        "lambda_z",
        "number",
        "shift",
        "nbproc",
        "w_n",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!OO!", kwlist, &PyArray_Type, &w_c_tmp, &PyArray_Type, &w_ghost_tmp, &PyArray_Type, &w_halo_tmp, &PyArray_Type, &w_haloghost_tmp, &PyArray_Type, &centerc_tmp, &PyArray_Type, &centerh_tmp, &PyArray_Type, &cellidn_tmp, &PyArray_Type, &periodicn_tmp, &PyArray_Type, &haloidn_tmp, &PyArray_Type, &vertexn_tmp, &PyArray_Type, &namen_tmp, &PyArray_Type, &centergn_tmp, &PyArray_Type, &halocentergn_tmp, &PyArray_Type, &R_x_tmp, &PyArray_Type, &R_y_tmp, &PyArray_Type, &R_z_tmp, &PyArray_Type, &lambda_x_tmp, &PyArray_Type, &lambda_y_tmp, &PyArray_Type, &lambda_z_tmp, &PyArray_Type, &number_tmp, &PyArray_Type, &shift_tmp, &nbproc_tmp, &PyArray_Type, &w_n_tmp))
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
    if (!pyarray_check(w_ghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_ghost = pyarray_to_ndarray(w_ghost_tmp);
    }
    if (!pyarray_check(w_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_halo = pyarray_to_ndarray(w_halo_tmp);
    }
    if (!pyarray_check(w_haloghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_haloghost = pyarray_to_ndarray(w_haloghost_tmp);
    }
    if (!pyarray_check(centerc_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerc = pyarray_to_ndarray(centerc_tmp);
    }
    if (!pyarray_check(centerh_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerh = pyarray_to_ndarray(centerh_tmp);
    }
    if (!pyarray_check(cellidn_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellidn = pyarray_to_ndarray(cellidn_tmp);
    }
    if (!pyarray_check(periodicn_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        periodicn = pyarray_to_ndarray(periodicn_tmp);
    }
    if (!pyarray_check(haloidn_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        haloidn = pyarray_to_ndarray(haloidn_tmp);
    }
    if (!pyarray_check(vertexn_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        vertexn = pyarray_to_ndarray(vertexn_tmp);
    }
    if (!pyarray_check(namen_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        namen = pyarray_to_ndarray(namen_tmp);
    }
    if (!pyarray_check(centergn_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centergn = pyarray_to_ndarray(centergn_tmp);
    }
    if (!pyarray_check(halocentergn_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        halocentergn = pyarray_to_ndarray(halocentergn_tmp);
    }
    if (!pyarray_check(R_x_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        R_x = pyarray_to_ndarray(R_x_tmp);
    }
    if (!pyarray_check(R_y_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        R_y = pyarray_to_ndarray(R_y_tmp);
    }
    if (!pyarray_check(R_z_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        R_z = pyarray_to_ndarray(R_z_tmp);
    }
    if (!pyarray_check(lambda_x_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        lambda_x = pyarray_to_ndarray(lambda_x_tmp);
    }
    if (!pyarray_check(lambda_y_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        lambda_y = pyarray_to_ndarray(lambda_y_tmp);
    }
    if (!pyarray_check(lambda_z_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        lambda_z = pyarray_to_ndarray(lambda_z_tmp);
    }
    if (!pyarray_check(number_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        number = pyarray_to_ndarray(number_tmp);
    }
    if (!pyarray_check(shift_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        shift = pyarray_to_ndarray(shift_tmp);
    }
    if (PyIs_NativeInt(nbproc_tmp))
    {
        nbproc = PyInt64_to_Int64(nbproc_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    if (!pyarray_check(w_n_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_n = pyarray_to_ndarray(w_n_tmp);
    }
    bind_c_centertovertex_3d(nd_ndim(&w_c, 0), nd_data(&w_c), nd_ndim(&w_ghost, 0), nd_data(&w_ghost), nd_ndim(&w_halo, 0), nd_data(&w_halo), nd_ndim(&w_haloghost, 0), nd_data(&w_haloghost), nd_ndim(&centerc, 0), nd_ndim(&centerc, 1), nd_data(&centerc), nd_ndim(&centerh, 0), nd_ndim(&centerh, 1), nd_data(&centerh), nd_ndim(&cellidn, 0), nd_ndim(&cellidn, 1), nd_data(&cellidn), nd_ndim(&periodicn, 0), nd_ndim(&periodicn, 1), nd_data(&periodicn), nd_ndim(&haloidn, 0), nd_ndim(&haloidn, 1), nd_data(&haloidn), nd_ndim(&vertexn, 0), nd_ndim(&vertexn, 1), nd_data(&vertexn), nd_ndim(&namen, 0), nd_data(&namen), nd_ndim(&centergn, 0), nd_ndim(&centergn, 1), nd_ndim(&centergn, 2), nd_data(&centergn), nd_ndim(&halocentergn, 0), nd_ndim(&halocentergn, 1), nd_ndim(&halocentergn, 2), nd_data(&halocentergn), nd_ndim(&R_x, 0), nd_data(&R_x), nd_ndim(&R_y, 0), nd_data(&R_y), nd_ndim(&R_z, 0), nd_data(&R_z), nd_ndim(&lambda_x, 0), nd_data(&lambda_x), nd_ndim(&lambda_y, 0), nd_data(&lambda_y), nd_ndim(&lambda_z, 0), nd_data(&lambda_z), nd_ndim(&number, 0), nd_data(&number), nd_ndim(&shift, 0), nd_ndim(&shift, 1), nd_data(&shift), nbproc, nd_ndim(&w_n, 0), nd_data(&w_n));
    result = Py_BuildValue("");
    free_pointer(w_c);
    free_pointer(w_ghost);
    free_pointer(w_halo);
    free_pointer(w_haloghost);
    free_pointer(centerc);
    free_pointer(centerh);
    free_pointer(cellidn);
    free_pointer(periodicn);
    free_pointer(haloidn);
    free_pointer(vertexn);
    free_pointer(namen);
    free_pointer(centergn);
    free_pointer(halocentergn);
    free_pointer(R_x);
    free_pointer(R_y);
    free_pointer(R_z);
    free_pointer(lambda_x);
    free_pointer(lambda_y);
    free_pointer(lambda_z);
    free_pointer(number);
    free_pointer(shift);
    free_pointer(w_n);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *barthlimiter_2d_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray w_c = {.shape = NULL};
    t_ndarray w_ghost = {.shape = NULL};
    t_ndarray w_halo = {.shape = NULL};
    t_ndarray w_x = {.shape = NULL};
    t_ndarray w_y = {.shape = NULL};
    t_ndarray w_z = {.shape = NULL};
    t_ndarray psi = {.shape = NULL};
    t_ndarray cellid = {.shape = NULL};
    t_ndarray faceid = {.shape = NULL};
    t_ndarray namef = {.shape = NULL};
    t_ndarray halofid = {.shape = NULL};
    t_ndarray centerc = {.shape = NULL};
    t_ndarray centerf = {.shape = NULL};
    PyArrayObject *w_c_tmp;
    PyArrayObject *w_ghost_tmp;
    PyArrayObject *w_halo_tmp;
    PyArrayObject *w_x_tmp;
    PyArrayObject *w_y_tmp;
    PyArrayObject *w_z_tmp;
    PyArrayObject *psi_tmp;
    PyArrayObject *cellid_tmp;
    PyArrayObject *faceid_tmp;
    PyArrayObject *namef_tmp;
    PyArrayObject *halofid_tmp;
    PyArrayObject *centerc_tmp;
    PyArrayObject *centerf_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "w_c",
        "w_ghost",
        "w_halo",
        "w_x",
        "w_y",
        "w_z",
        "psi",
        "cellid",
        "faceid",
        "namef",
        "halofid",
        "centerc",
        "centerf",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!O!O!O!O!O!O!", kwlist, &PyArray_Type, &w_c_tmp, &PyArray_Type, &w_ghost_tmp, &PyArray_Type, &w_halo_tmp, &PyArray_Type, &w_x_tmp, &PyArray_Type, &w_y_tmp, &PyArray_Type, &w_z_tmp, &PyArray_Type, &psi_tmp, &PyArray_Type, &cellid_tmp, &PyArray_Type, &faceid_tmp, &PyArray_Type, &namef_tmp, &PyArray_Type, &halofid_tmp, &PyArray_Type, &centerc_tmp, &PyArray_Type, &centerf_tmp))
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
    if (!pyarray_check(w_ghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_ghost = pyarray_to_ndarray(w_ghost_tmp);
    }
    if (!pyarray_check(w_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_halo = pyarray_to_ndarray(w_halo_tmp);
    }
    if (!pyarray_check(w_x_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_x = pyarray_to_ndarray(w_x_tmp);
    }
    if (!pyarray_check(w_y_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_y = pyarray_to_ndarray(w_y_tmp);
    }
    if (!pyarray_check(w_z_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        w_z = pyarray_to_ndarray(w_z_tmp);
    }
    if (!pyarray_check(psi_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        psi = pyarray_to_ndarray(psi_tmp);
    }
    if (!pyarray_check(cellid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellid = pyarray_to_ndarray(cellid_tmp);
    }
    if (!pyarray_check(faceid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        faceid = pyarray_to_ndarray(faceid_tmp);
    }
    if (!pyarray_check(namef_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        namef = pyarray_to_ndarray(namef_tmp);
    }
    if (!pyarray_check(halofid_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofid = pyarray_to_ndarray(halofid_tmp);
    }
    if (!pyarray_check(centerc_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerc = pyarray_to_ndarray(centerc_tmp);
    }
    if (!pyarray_check(centerf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerf = pyarray_to_ndarray(centerf_tmp);
    }
    bind_c_barthlimiter_2d(nd_ndim(&w_c, 0), nd_data(&w_c), nd_ndim(&w_ghost, 0), nd_data(&w_ghost), nd_ndim(&w_halo, 0), nd_data(&w_halo), nd_ndim(&w_x, 0), nd_data(&w_x), nd_ndim(&w_y, 0), nd_data(&w_y), nd_ndim(&w_z, 0), nd_data(&w_z), nd_ndim(&psi, 0), nd_data(&psi), nd_ndim(&cellid, 0), nd_ndim(&cellid, 1), nd_data(&cellid), nd_ndim(&faceid, 0), nd_ndim(&faceid, 1), nd_data(&faceid), nd_ndim(&namef, 0), nd_data(&namef), nd_ndim(&halofid, 0), nd_data(&halofid), nd_ndim(&centerc, 0), nd_ndim(&centerc, 1), nd_data(&centerc), nd_ndim(&centerf, 0), nd_ndim(&centerf, 1), nd_data(&centerf));
    result = Py_BuildValue("");
    free_pointer(w_c);
    free_pointer(w_ghost);
    free_pointer(w_halo);
    free_pointer(w_x);
    free_pointer(w_y);
    free_pointer(w_z);
    free_pointer(psi);
    free_pointer(cellid);
    free_pointer(faceid);
    free_pointer(namef);
    free_pointer(halofid);
    free_pointer(centerc);
    free_pointer(centerf);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *barthlimiter_3d_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray h_c = {.shape = NULL};
    t_ndarray h_ghost = {.shape = NULL};
    t_ndarray h_halo = {.shape = NULL};
    t_ndarray h_x = {.shape = NULL};
    t_ndarray h_y = {.shape = NULL};
    t_ndarray h_z = {.shape = NULL};
    t_ndarray psi = {.shape = NULL};
    t_ndarray cellid = {.shape = NULL};
    t_ndarray faceid = {.shape = NULL};
    t_ndarray namef = {.shape = NULL};
    t_ndarray halofid = {.shape = NULL};
    t_ndarray centerc = {.shape = NULL};
    t_ndarray centerf = {.shape = NULL};
    PyArrayObject *h_c_tmp;
    PyArrayObject *h_ghost_tmp;
    PyArrayObject *h_halo_tmp;
    PyArrayObject *h_x_tmp;
    PyArrayObject *h_y_tmp;
    PyArrayObject *h_z_tmp;
    PyArrayObject *psi_tmp;
    PyArrayObject *cellid_tmp;
    PyArrayObject *faceid_tmp;
    PyArrayObject *namef_tmp;
    PyArrayObject *halofid_tmp;
    PyArrayObject *centerc_tmp;
    PyArrayObject *centerf_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "h_c",
        "h_ghost",
        "h_halo",
        "h_x",
        "h_y",
        "h_z",
        "psi",
        "cellid",
        "faceid",
        "namef",
        "halofid",
        "centerc",
        "centerf",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!O!O!O!O!O!O!", kwlist, &PyArray_Type, &h_c_tmp, &PyArray_Type, &h_ghost_tmp, &PyArray_Type, &h_halo_tmp, &PyArray_Type, &h_x_tmp, &PyArray_Type, &h_y_tmp, &PyArray_Type, &h_z_tmp, &PyArray_Type, &psi_tmp, &PyArray_Type, &cellid_tmp, &PyArray_Type, &faceid_tmp, &PyArray_Type, &namef_tmp, &PyArray_Type, &halofid_tmp, &PyArray_Type, &centerc_tmp, &PyArray_Type, &centerf_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(h_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        h_c = pyarray_to_ndarray(h_c_tmp);
    }
    if (!pyarray_check(h_ghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        h_ghost = pyarray_to_ndarray(h_ghost_tmp);
    }
    if (!pyarray_check(h_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        h_halo = pyarray_to_ndarray(h_halo_tmp);
    }
    if (!pyarray_check(h_x_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        h_x = pyarray_to_ndarray(h_x_tmp);
    }
    if (!pyarray_check(h_y_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        h_y = pyarray_to_ndarray(h_y_tmp);
    }
    if (!pyarray_check(h_z_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        h_z = pyarray_to_ndarray(h_z_tmp);
    }
    if (!pyarray_check(psi_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        psi = pyarray_to_ndarray(psi_tmp);
    }
    if (!pyarray_check(cellid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellid = pyarray_to_ndarray(cellid_tmp);
    }
    if (!pyarray_check(faceid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        faceid = pyarray_to_ndarray(faceid_tmp);
    }
    if (!pyarray_check(namef_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        namef = pyarray_to_ndarray(namef_tmp);
    }
    if (!pyarray_check(halofid_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofid = pyarray_to_ndarray(halofid_tmp);
    }
    if (!pyarray_check(centerc_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerc = pyarray_to_ndarray(centerc_tmp);
    }
    if (!pyarray_check(centerf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerf = pyarray_to_ndarray(centerf_tmp);
    }
    bind_c_barthlimiter_3d(nd_ndim(&h_c, 0), nd_data(&h_c), nd_ndim(&h_ghost, 0), nd_data(&h_ghost), nd_ndim(&h_halo, 0), nd_data(&h_halo), nd_ndim(&h_x, 0), nd_data(&h_x), nd_ndim(&h_y, 0), nd_data(&h_y), nd_ndim(&h_z, 0), nd_data(&h_z), nd_ndim(&psi, 0), nd_data(&psi), nd_ndim(&cellid, 0), nd_ndim(&cellid, 1), nd_data(&cellid), nd_ndim(&faceid, 0), nd_ndim(&faceid, 1), nd_data(&faceid), nd_ndim(&namef, 0), nd_data(&namef), nd_ndim(&halofid, 0), nd_data(&halofid), nd_ndim(&centerc, 0), nd_ndim(&centerc, 1), nd_data(&centerc), nd_ndim(&centerf, 0), nd_ndim(&centerf, 1), nd_data(&centerf));
    result = Py_BuildValue("");
    free_pointer(h_c);
    free_pointer(h_ghost);
    free_pointer(h_halo);
    free_pointer(h_x);
    free_pointer(h_y);
    free_pointer(h_z);
    free_pointer(psi);
    free_pointer(cellid);
    free_pointer(faceid);
    free_pointer(namef);
    free_pointer(halofid);
    free_pointer(centerc);
    free_pointer(centerf);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *search_element_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray a = {.shape = NULL};
    int64_t target_value;
    int64_t find;
    PyArrayObject *a_tmp;
    PyObject *target_value_tmp;
    PyObject *find_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "a",
        "target_value",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O", kwlist, &PyArray_Type, &a_tmp, &target_value_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(a_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        a = pyarray_to_ndarray(a_tmp);
    }
    if (PyIs_NativeInt(target_value_tmp))
    {
        target_value = PyInt64_to_Int64(target_value_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    find = bind_c_search_element(nd_ndim(&a, 0), nd_data(&a), target_value);
    find_tmp = Int64_to_PyLong(&find);
    result = Py_BuildValue("O", find_tmp);
    Py_DECREF(find_tmp);
    free_pointer(a);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *get_triplet_2d_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray cellfid = {.shape = NULL};
    t_ndarray nodeidf = {.shape = NULL};
    t_ndarray vertexn = {.shape = NULL};
    t_ndarray halofid = {.shape = NULL};
    t_ndarray haloext = {.shape = NULL};
    t_ndarray namen = {.shape = NULL};
    t_ndarray oldnamen = {.shape = NULL};
    t_ndarray volume = {.shape = NULL};
    t_ndarray cellnid = {.shape = NULL};
    t_ndarray centerc = {.shape = NULL};
    t_ndarray centerh = {.shape = NULL};
    t_ndarray halonid = {.shape = NULL};
    t_ndarray periodicnid = {.shape = NULL};
    t_ndarray centergn = {.shape = NULL};
    t_ndarray halocentergn = {.shape = NULL};
    t_ndarray airDiamond = {.shape = NULL};
    t_ndarray lambda_x = {.shape = NULL};
    t_ndarray lambda_y = {.shape = NULL};
    t_ndarray number = {.shape = NULL};
    t_ndarray R_x = {.shape = NULL};
    t_ndarray R_y = {.shape = NULL};
    t_ndarray param1 = {.shape = NULL};
    t_ndarray param2 = {.shape = NULL};
    t_ndarray param3 = {.shape = NULL};
    t_ndarray param4 = {.shape = NULL};
    t_ndarray shift = {.shape = NULL};
    int64_t nbelements;
    t_ndarray loctoglob = {.shape = NULL};
    t_ndarray BCdirichlet = {.shape = NULL};
    t_ndarray a_loc = {.shape = NULL};
    t_ndarray irn_loc = {.shape = NULL};
    t_ndarray jcn_loc = {.shape = NULL};
    t_ndarray matrixinnerfaces = {.shape = NULL};
    t_ndarray halofaces = {.shape = NULL};
    t_ndarray dirichletfaces = {.shape = NULL};
    PyArrayObject *cellfid_tmp;
    PyArrayObject *nodeidf_tmp;
    PyArrayObject *vertexn_tmp;
    PyArrayObject *halofid_tmp;
    PyArrayObject *haloext_tmp;
    PyArrayObject *namen_tmp;
    PyArrayObject *oldnamen_tmp;
    PyArrayObject *volume_tmp;
    PyArrayObject *cellnid_tmp;
    PyArrayObject *centerc_tmp;
    PyArrayObject *centerh_tmp;
    PyArrayObject *halonid_tmp;
    PyArrayObject *periodicnid_tmp;
    PyArrayObject *centergn_tmp;
    PyArrayObject *halocentergn_tmp;
    PyArrayObject *airDiamond_tmp;
    PyArrayObject *lambda_x_tmp;
    PyArrayObject *lambda_y_tmp;
    PyArrayObject *number_tmp;
    PyArrayObject *R_x_tmp;
    PyArrayObject *R_y_tmp;
    PyArrayObject *param1_tmp;
    PyArrayObject *param2_tmp;
    PyArrayObject *param3_tmp;
    PyArrayObject *param4_tmp;
    PyArrayObject *shift_tmp;
    PyObject *nbelements_tmp;
    PyArrayObject *loctoglob_tmp;
    PyArrayObject *BCdirichlet_tmp;
    PyArrayObject *a_loc_tmp;
    PyArrayObject *irn_loc_tmp;
    PyArrayObject *jcn_loc_tmp;
    PyArrayObject *matrixinnerfaces_tmp;
    PyArrayObject *halofaces_tmp;
    PyArrayObject *dirichletfaces_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "cellfid",
        "nodeidf",
        "vertexn",
        "halofid",
        "haloext",
        "namen",
        "oldnamen",
        "volume",
        "cellnid",
        "centerc",
        "centerh",
        "halonid",
        "periodicnid",
        "centergn",
        "halocentergn",
        "airDiamond",
        "lambda_x",
        "lambda_y",
        "number",
        "R_x",
        "R_y",
        "param1",
        "param2",
        "param3",
        "param4",
        "shift",
        "nbelements",
        "loctoglob",
        "BCdirichlet",
        "a_loc",
        "irn_loc",
        "jcn_loc",
        "matrixinnerfaces",
        "halofaces",
        "dirichletfaces",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!OO!O!O!O!O!O!O!O!", kwlist, &PyArray_Type, &cellfid_tmp, &PyArray_Type, &nodeidf_tmp, &PyArray_Type, &vertexn_tmp, &PyArray_Type, &halofid_tmp, &PyArray_Type, &haloext_tmp, &PyArray_Type, &namen_tmp, &PyArray_Type, &oldnamen_tmp, &PyArray_Type, &volume_tmp, &PyArray_Type, &cellnid_tmp, &PyArray_Type, &centerc_tmp, &PyArray_Type, &centerh_tmp, &PyArray_Type, &halonid_tmp, &PyArray_Type, &periodicnid_tmp, &PyArray_Type, &centergn_tmp, &PyArray_Type, &halocentergn_tmp, &PyArray_Type, &airDiamond_tmp, &PyArray_Type, &lambda_x_tmp, &PyArray_Type, &lambda_y_tmp, &PyArray_Type, &number_tmp, &PyArray_Type, &R_x_tmp, &PyArray_Type, &R_y_tmp, &PyArray_Type, &param1_tmp, &PyArray_Type, &param2_tmp, &PyArray_Type, &param3_tmp, &PyArray_Type, &param4_tmp, &PyArray_Type, &shift_tmp, &nbelements_tmp, &PyArray_Type, &loctoglob_tmp, &PyArray_Type, &BCdirichlet_tmp, &PyArray_Type, &a_loc_tmp, &PyArray_Type, &irn_loc_tmp, &PyArray_Type, &jcn_loc_tmp, &PyArray_Type, &matrixinnerfaces_tmp, &PyArray_Type, &halofaces_tmp, &PyArray_Type, &dirichletfaces_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(cellfid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellfid = pyarray_to_ndarray(cellfid_tmp);
    }
    if (!pyarray_check(nodeidf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nodeidf = pyarray_to_ndarray(nodeidf_tmp);
    }
    if (!pyarray_check(vertexn_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        vertexn = pyarray_to_ndarray(vertexn_tmp);
    }
    if (!pyarray_check(halofid_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofid = pyarray_to_ndarray(halofid_tmp);
    }
    if (!pyarray_check(haloext_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        haloext = pyarray_to_ndarray(haloext_tmp);
    }
    if (!pyarray_check(namen_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        namen = pyarray_to_ndarray(namen_tmp);
    }
    if (!pyarray_check(oldnamen_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        oldnamen = pyarray_to_ndarray(oldnamen_tmp);
    }
    if (!pyarray_check(volume_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        volume = pyarray_to_ndarray(volume_tmp);
    }
    if (!pyarray_check(cellnid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellnid = pyarray_to_ndarray(cellnid_tmp);
    }
    if (!pyarray_check(centerc_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerc = pyarray_to_ndarray(centerc_tmp);
    }
    if (!pyarray_check(centerh_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerh = pyarray_to_ndarray(centerh_tmp);
    }
    if (!pyarray_check(halonid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        halonid = pyarray_to_ndarray(halonid_tmp);
    }
    if (!pyarray_check(periodicnid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        periodicnid = pyarray_to_ndarray(periodicnid_tmp);
    }
    if (!pyarray_check(centergn_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centergn = pyarray_to_ndarray(centergn_tmp);
    }
    if (!pyarray_check(halocentergn_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        halocentergn = pyarray_to_ndarray(halocentergn_tmp);
    }
    if (!pyarray_check(airDiamond_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        airDiamond = pyarray_to_ndarray(airDiamond_tmp);
    }
    if (!pyarray_check(lambda_x_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        lambda_x = pyarray_to_ndarray(lambda_x_tmp);
    }
    if (!pyarray_check(lambda_y_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        lambda_y = pyarray_to_ndarray(lambda_y_tmp);
    }
    if (!pyarray_check(number_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        number = pyarray_to_ndarray(number_tmp);
    }
    if (!pyarray_check(R_x_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        R_x = pyarray_to_ndarray(R_x_tmp);
    }
    if (!pyarray_check(R_y_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        R_y = pyarray_to_ndarray(R_y_tmp);
    }
    if (!pyarray_check(param1_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        param1 = pyarray_to_ndarray(param1_tmp);
    }
    if (!pyarray_check(param2_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        param2 = pyarray_to_ndarray(param2_tmp);
    }
    if (!pyarray_check(param3_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        param3 = pyarray_to_ndarray(param3_tmp);
    }
    if (!pyarray_check(param4_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        param4 = pyarray_to_ndarray(param4_tmp);
    }
    if (!pyarray_check(shift_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        shift = pyarray_to_ndarray(shift_tmp);
    }
    if (PyIs_NativeInt(nbelements_tmp))
    {
        nbelements = PyInt64_to_Int64(nbelements_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    if (!pyarray_check(loctoglob_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        loctoglob = pyarray_to_ndarray(loctoglob_tmp);
    }
    if (!pyarray_check(BCdirichlet_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        BCdirichlet = pyarray_to_ndarray(BCdirichlet_tmp);
    }
    if (!pyarray_check(a_loc_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        a_loc = pyarray_to_ndarray(a_loc_tmp);
    }
    if (!pyarray_check(irn_loc_tmp, NPY_INT32, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        irn_loc = pyarray_to_ndarray(irn_loc_tmp);
    }
    if (!pyarray_check(jcn_loc_tmp, NPY_INT32, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        jcn_loc = pyarray_to_ndarray(jcn_loc_tmp);
    }
    if (!pyarray_check(matrixinnerfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        matrixinnerfaces = pyarray_to_ndarray(matrixinnerfaces_tmp);
    }
    if (!pyarray_check(halofaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofaces = pyarray_to_ndarray(halofaces_tmp);
    }
    if (!pyarray_check(dirichletfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        dirichletfaces = pyarray_to_ndarray(dirichletfaces_tmp);
    }
    bind_c_get_triplet_2d(nd_ndim(&cellfid, 0), nd_ndim(&cellfid, 1), nd_data(&cellfid), nd_ndim(&nodeidf, 0), nd_ndim(&nodeidf, 1), nd_data(&nodeidf), nd_ndim(&vertexn, 0), nd_ndim(&vertexn, 1), nd_data(&vertexn), nd_ndim(&halofid, 0), nd_data(&halofid), nd_ndim(&haloext, 0), nd_ndim(&haloext, 1), nd_data(&haloext), nd_ndim(&namen, 0), nd_data(&namen), nd_ndim(&oldnamen, 0), nd_data(&oldnamen), nd_ndim(&volume, 0), nd_data(&volume), nd_ndim(&cellnid, 0), nd_ndim(&cellnid, 1), nd_data(&cellnid), nd_ndim(&centerc, 0), nd_ndim(&centerc, 1), nd_data(&centerc), nd_ndim(&centerh, 0), nd_ndim(&centerh, 1), nd_data(&centerh), nd_ndim(&halonid, 0), nd_ndim(&halonid, 1), nd_data(&halonid), nd_ndim(&periodicnid, 0), nd_ndim(&periodicnid, 1), nd_data(&periodicnid), nd_ndim(&centergn, 0), nd_ndim(&centergn, 1), nd_ndim(&centergn, 2), nd_data(&centergn), nd_ndim(&halocentergn, 0), nd_ndim(&halocentergn, 1), nd_ndim(&halocentergn, 2), nd_data(&halocentergn), nd_ndim(&airDiamond, 0), nd_data(&airDiamond), nd_ndim(&lambda_x, 0), nd_data(&lambda_x), nd_ndim(&lambda_y, 0), nd_data(&lambda_y), nd_ndim(&number, 0), nd_data(&number), nd_ndim(&R_x, 0), nd_data(&R_x), nd_ndim(&R_y, 0), nd_data(&R_y), nd_ndim(&param1, 0), nd_data(&param1), nd_ndim(&param2, 0), nd_data(&param2), nd_ndim(&param3, 0), nd_data(&param3), nd_ndim(&param4, 0), nd_data(&param4), nd_ndim(&shift, 0), nd_ndim(&shift, 1), nd_data(&shift), nbelements, nd_ndim(&loctoglob, 0), nd_data(&loctoglob), nd_ndim(&BCdirichlet, 0), nd_data(&BCdirichlet), nd_ndim(&a_loc, 0), nd_data(&a_loc), nd_ndim(&irn_loc, 0), nd_data(&irn_loc), nd_ndim(&jcn_loc, 0), nd_data(&jcn_loc), nd_ndim(&matrixinnerfaces, 0), nd_data(&matrixinnerfaces), nd_ndim(&halofaces, 0), nd_data(&halofaces), nd_ndim(&dirichletfaces, 0), nd_data(&dirichletfaces));
    result = Py_BuildValue("");
    free_pointer(cellfid);
    free_pointer(nodeidf);
    free_pointer(vertexn);
    free_pointer(halofid);
    free_pointer(haloext);
    free_pointer(namen);
    free_pointer(oldnamen);
    free_pointer(volume);
    free_pointer(cellnid);
    free_pointer(centerc);
    free_pointer(centerh);
    free_pointer(halonid);
    free_pointer(periodicnid);
    free_pointer(centergn);
    free_pointer(halocentergn);
    free_pointer(airDiamond);
    free_pointer(lambda_x);
    free_pointer(lambda_y);
    free_pointer(number);
    free_pointer(R_x);
    free_pointer(R_y);
    free_pointer(param1);
    free_pointer(param2);
    free_pointer(param3);
    free_pointer(param4);
    free_pointer(shift);
    free_pointer(loctoglob);
    free_pointer(BCdirichlet);
    free_pointer(a_loc);
    free_pointer(irn_loc);
    free_pointer(jcn_loc);
    free_pointer(matrixinnerfaces);
    free_pointer(halofaces);
    free_pointer(dirichletfaces);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *compute_2dmatrix_size_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray nodeidf = {.shape = NULL};
    t_ndarray halofid = {.shape = NULL};
    t_ndarray cellnid = {.shape = NULL};
    t_ndarray halonid = {.shape = NULL};
    t_ndarray periodicnid = {.shape = NULL};
    t_ndarray centergn = {.shape = NULL};
    t_ndarray halocentergn = {.shape = NULL};
    t_ndarray oldnamen = {.shape = NULL};
    t_ndarray BCdirichlet = {.shape = NULL};
    t_ndarray matrixinnerfaces = {.shape = NULL};
    t_ndarray halofaces = {.shape = NULL};
    t_ndarray dirichletfaces = {.shape = NULL};
    int64_t cmpt;
    PyArrayObject *nodeidf_tmp;
    PyArrayObject *halofid_tmp;
    PyArrayObject *cellnid_tmp;
    PyArrayObject *halonid_tmp;
    PyArrayObject *periodicnid_tmp;
    PyArrayObject *centergn_tmp;
    PyArrayObject *halocentergn_tmp;
    PyArrayObject *oldnamen_tmp;
    PyArrayObject *BCdirichlet_tmp;
    PyArrayObject *matrixinnerfaces_tmp;
    PyArrayObject *halofaces_tmp;
    PyArrayObject *dirichletfaces_tmp;
    PyObject *cmpt_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "nodeidf",
        "halofid",
        "cellnid",
        "halonid",
        "periodicnid",
        "centergn",
        "halocentergn",
        "oldnamen",
        "BCdirichlet",
        "matrixinnerfaces",
        "halofaces",
        "dirichletfaces",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!O!O!O!O!O!", kwlist, &PyArray_Type, &nodeidf_tmp, &PyArray_Type, &halofid_tmp, &PyArray_Type, &cellnid_tmp, &PyArray_Type, &halonid_tmp, &PyArray_Type, &periodicnid_tmp, &PyArray_Type, &centergn_tmp, &PyArray_Type, &halocentergn_tmp, &PyArray_Type, &oldnamen_tmp, &PyArray_Type, &BCdirichlet_tmp, &PyArray_Type, &matrixinnerfaces_tmp, &PyArray_Type, &halofaces_tmp, &PyArray_Type, &dirichletfaces_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(nodeidf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nodeidf = pyarray_to_ndarray(nodeidf_tmp);
    }
    if (!pyarray_check(halofid_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofid = pyarray_to_ndarray(halofid_tmp);
    }
    if (!pyarray_check(cellnid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellnid = pyarray_to_ndarray(cellnid_tmp);
    }
    if (!pyarray_check(halonid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        halonid = pyarray_to_ndarray(halonid_tmp);
    }
    if (!pyarray_check(periodicnid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        periodicnid = pyarray_to_ndarray(periodicnid_tmp);
    }
    if (!pyarray_check(centergn_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centergn = pyarray_to_ndarray(centergn_tmp);
    }
    if (!pyarray_check(halocentergn_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        halocentergn = pyarray_to_ndarray(halocentergn_tmp);
    }
    if (!pyarray_check(oldnamen_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        oldnamen = pyarray_to_ndarray(oldnamen_tmp);
    }
    if (!pyarray_check(BCdirichlet_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        BCdirichlet = pyarray_to_ndarray(BCdirichlet_tmp);
    }
    if (!pyarray_check(matrixinnerfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        matrixinnerfaces = pyarray_to_ndarray(matrixinnerfaces_tmp);
    }
    if (!pyarray_check(halofaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofaces = pyarray_to_ndarray(halofaces_tmp);
    }
    if (!pyarray_check(dirichletfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        dirichletfaces = pyarray_to_ndarray(dirichletfaces_tmp);
    }
    cmpt = bind_c_compute_2dmatrix_size(nd_ndim(&nodeidf, 0), nd_ndim(&nodeidf, 1), nd_data(&nodeidf), nd_ndim(&halofid, 0), nd_data(&halofid), nd_ndim(&cellnid, 0), nd_ndim(&cellnid, 1), nd_data(&cellnid), nd_ndim(&halonid, 0), nd_ndim(&halonid, 1), nd_data(&halonid), nd_ndim(&periodicnid, 0), nd_ndim(&periodicnid, 1), nd_data(&periodicnid), nd_ndim(&centergn, 0), nd_ndim(&centergn, 1), nd_ndim(&centergn, 2), nd_data(&centergn), nd_ndim(&halocentergn, 0), nd_ndim(&halocentergn, 1), nd_ndim(&halocentergn, 2), nd_data(&halocentergn), nd_ndim(&oldnamen, 0), nd_data(&oldnamen), nd_ndim(&BCdirichlet, 0), nd_data(&BCdirichlet), nd_ndim(&matrixinnerfaces, 0), nd_data(&matrixinnerfaces), nd_ndim(&halofaces, 0), nd_data(&halofaces), nd_ndim(&dirichletfaces, 0), nd_data(&dirichletfaces));
    cmpt_tmp = Int64_to_PyLong(&cmpt);
    result = Py_BuildValue("O", cmpt_tmp);
    Py_DECREF(cmpt_tmp);
    free_pointer(nodeidf);
    free_pointer(halofid);
    free_pointer(cellnid);
    free_pointer(halonid);
    free_pointer(periodicnid);
    free_pointer(centergn);
    free_pointer(halocentergn);
    free_pointer(oldnamen);
    free_pointer(BCdirichlet);
    free_pointer(matrixinnerfaces);
    free_pointer(halofaces);
    free_pointer(dirichletfaces);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *compute_3dmatrix_size_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray nodeidf = {.shape = NULL};
    t_ndarray halofid = {.shape = NULL};
    t_ndarray cellnid = {.shape = NULL};
    t_ndarray halonid = {.shape = NULL};
    t_ndarray periodicnid = {.shape = NULL};
    t_ndarray centergn = {.shape = NULL};
    t_ndarray halocentergn = {.shape = NULL};
    t_ndarray oldnamen = {.shape = NULL};
    t_ndarray BCdirichlet = {.shape = NULL};
    t_ndarray matrixinnerfaces = {.shape = NULL};
    t_ndarray halofaces = {.shape = NULL};
    t_ndarray dirichletfaces = {.shape = NULL};
    int64_t cmpt;
    PyArrayObject *nodeidf_tmp;
    PyArrayObject *halofid_tmp;
    PyArrayObject *cellnid_tmp;
    PyArrayObject *halonid_tmp;
    PyArrayObject *periodicnid_tmp;
    PyArrayObject *centergn_tmp;
    PyArrayObject *halocentergn_tmp;
    PyArrayObject *oldnamen_tmp;
    PyArrayObject *BCdirichlet_tmp;
    PyArrayObject *matrixinnerfaces_tmp;
    PyArrayObject *halofaces_tmp;
    PyArrayObject *dirichletfaces_tmp;
    PyObject *cmpt_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "nodeidf",
        "halofid",
        "cellnid",
        "halonid",
        "periodicnid",
        "centergn",
        "halocentergn",
        "oldnamen",
        "BCdirichlet",
        "matrixinnerfaces",
        "halofaces",
        "dirichletfaces",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!O!O!O!O!O!", kwlist, &PyArray_Type, &nodeidf_tmp, &PyArray_Type, &halofid_tmp, &PyArray_Type, &cellnid_tmp, &PyArray_Type, &halonid_tmp, &PyArray_Type, &periodicnid_tmp, &PyArray_Type, &centergn_tmp, &PyArray_Type, &halocentergn_tmp, &PyArray_Type, &oldnamen_tmp, &PyArray_Type, &BCdirichlet_tmp, &PyArray_Type, &matrixinnerfaces_tmp, &PyArray_Type, &halofaces_tmp, &PyArray_Type, &dirichletfaces_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(nodeidf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nodeidf = pyarray_to_ndarray(nodeidf_tmp);
    }
    if (!pyarray_check(halofid_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofid = pyarray_to_ndarray(halofid_tmp);
    }
    if (!pyarray_check(cellnid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellnid = pyarray_to_ndarray(cellnid_tmp);
    }
    if (!pyarray_check(halonid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        halonid = pyarray_to_ndarray(halonid_tmp);
    }
    if (!pyarray_check(periodicnid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        periodicnid = pyarray_to_ndarray(periodicnid_tmp);
    }
    if (!pyarray_check(centergn_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centergn = pyarray_to_ndarray(centergn_tmp);
    }
    if (!pyarray_check(halocentergn_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        halocentergn = pyarray_to_ndarray(halocentergn_tmp);
    }
    if (!pyarray_check(oldnamen_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        oldnamen = pyarray_to_ndarray(oldnamen_tmp);
    }
    if (!pyarray_check(BCdirichlet_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        BCdirichlet = pyarray_to_ndarray(BCdirichlet_tmp);
    }
    if (!pyarray_check(matrixinnerfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        matrixinnerfaces = pyarray_to_ndarray(matrixinnerfaces_tmp);
    }
    if (!pyarray_check(halofaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofaces = pyarray_to_ndarray(halofaces_tmp);
    }
    if (!pyarray_check(dirichletfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        dirichletfaces = pyarray_to_ndarray(dirichletfaces_tmp);
    }
    cmpt = bind_c_compute_3dmatrix_size(nd_ndim(&nodeidf, 0), nd_ndim(&nodeidf, 1), nd_data(&nodeidf), nd_ndim(&halofid, 0), nd_data(&halofid), nd_ndim(&cellnid, 0), nd_ndim(&cellnid, 1), nd_data(&cellnid), nd_ndim(&halonid, 0), nd_ndim(&halonid, 1), nd_data(&halonid), nd_ndim(&periodicnid, 0), nd_ndim(&periodicnid, 1), nd_data(&periodicnid), nd_ndim(&centergn, 0), nd_ndim(&centergn, 1), nd_ndim(&centergn, 2), nd_data(&centergn), nd_ndim(&halocentergn, 0), nd_ndim(&halocentergn, 1), nd_ndim(&halocentergn, 2), nd_data(&halocentergn), nd_ndim(&oldnamen, 0), nd_data(&oldnamen), nd_ndim(&BCdirichlet, 0), nd_data(&BCdirichlet), nd_ndim(&matrixinnerfaces, 0), nd_data(&matrixinnerfaces), nd_ndim(&halofaces, 0), nd_data(&halofaces), nd_ndim(&dirichletfaces, 0), nd_data(&dirichletfaces));
    cmpt_tmp = Int64_to_PyLong(&cmpt);
    result = Py_BuildValue("O", cmpt_tmp);
    Py_DECREF(cmpt_tmp);
    free_pointer(nodeidf);
    free_pointer(halofid);
    free_pointer(cellnid);
    free_pointer(halonid);
    free_pointer(periodicnid);
    free_pointer(centergn);
    free_pointer(halocentergn);
    free_pointer(oldnamen);
    free_pointer(BCdirichlet);
    free_pointer(matrixinnerfaces);
    free_pointer(halofaces);
    free_pointer(dirichletfaces);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *get_triplet_3d_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray cellfid = {.shape = NULL};
    t_ndarray nodeidf = {.shape = NULL};
    t_ndarray vertexn = {.shape = NULL};
    t_ndarray halofid = {.shape = NULL};
    t_ndarray haloext = {.shape = NULL};
    t_ndarray namen = {.shape = NULL};
    t_ndarray oldnamen = {.shape = NULL};
    t_ndarray volume = {.shape = NULL};
    t_ndarray centergn = {.shape = NULL};
    t_ndarray halocentergn = {.shape = NULL};
    t_ndarray periodicnid = {.shape = NULL};
    t_ndarray cellnid = {.shape = NULL};
    t_ndarray centerc = {.shape = NULL};
    t_ndarray centerh = {.shape = NULL};
    t_ndarray halonid = {.shape = NULL};
    t_ndarray airDiamond = {.shape = NULL};
    t_ndarray lambda_x = {.shape = NULL};
    t_ndarray lambda_y = {.shape = NULL};
    t_ndarray lambda_z = {.shape = NULL};
    t_ndarray number = {.shape = NULL};
    t_ndarray R_x = {.shape = NULL};
    t_ndarray R_y = {.shape = NULL};
    t_ndarray R_z = {.shape = NULL};
    t_ndarray param1 = {.shape = NULL};
    t_ndarray param2 = {.shape = NULL};
    t_ndarray param3 = {.shape = NULL};
    t_ndarray shift = {.shape = NULL};
    t_ndarray loctoglob = {.shape = NULL};
    t_ndarray BCdirichlet = {.shape = NULL};
    t_ndarray a_loc = {.shape = NULL};
    t_ndarray irn_loc = {.shape = NULL};
    t_ndarray jcn_loc = {.shape = NULL};
    t_ndarray matrixinnerfaces = {.shape = NULL};
    t_ndarray halofaces = {.shape = NULL};
    t_ndarray dirichletfaces = {.shape = NULL};
    PyArrayObject *cellfid_tmp;
    PyArrayObject *nodeidf_tmp;
    PyArrayObject *vertexn_tmp;
    PyArrayObject *halofid_tmp;
    PyArrayObject *haloext_tmp;
    PyArrayObject *namen_tmp;
    PyArrayObject *oldnamen_tmp;
    PyArrayObject *volume_tmp;
    PyArrayObject *centergn_tmp;
    PyArrayObject *halocentergn_tmp;
    PyArrayObject *periodicnid_tmp;
    PyArrayObject *cellnid_tmp;
    PyArrayObject *centerc_tmp;
    PyArrayObject *centerh_tmp;
    PyArrayObject *halonid_tmp;
    PyArrayObject *airDiamond_tmp;
    PyArrayObject *lambda_x_tmp;
    PyArrayObject *lambda_y_tmp;
    PyArrayObject *lambda_z_tmp;
    PyArrayObject *number_tmp;
    PyArrayObject *R_x_tmp;
    PyArrayObject *R_y_tmp;
    PyArrayObject *R_z_tmp;
    PyArrayObject *param1_tmp;
    PyArrayObject *param2_tmp;
    PyArrayObject *param3_tmp;
    PyArrayObject *shift_tmp;
    PyArrayObject *loctoglob_tmp;
    PyArrayObject *BCdirichlet_tmp;
    PyArrayObject *a_loc_tmp;
    PyArrayObject *irn_loc_tmp;
    PyArrayObject *jcn_loc_tmp;
    PyArrayObject *matrixinnerfaces_tmp;
    PyArrayObject *halofaces_tmp;
    PyArrayObject *dirichletfaces_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "cellfid",
        "nodeidf",
        "vertexn",
        "halofid",
        "haloext",
        "namen",
        "oldnamen",
        "volume",
        "centergn",
        "halocentergn",
        "periodicnid",
        "cellnid",
        "centerc",
        "centerh",
        "halonid",
        "airDiamond",
        "lambda_x",
        "lambda_y",
        "lambda_z",
        "number",
        "R_x",
        "R_y",
        "R_z",
        "param1",
        "param2",
        "param3",
        "shift",
        "loctoglob",
        "BCdirichlet",
        "a_loc",
        "irn_loc",
        "jcn_loc",
        "matrixinnerfaces",
        "halofaces",
        "dirichletfaces",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!", kwlist, &PyArray_Type, &cellfid_tmp, &PyArray_Type, &nodeidf_tmp, &PyArray_Type, &vertexn_tmp, &PyArray_Type, &halofid_tmp, &PyArray_Type, &haloext_tmp, &PyArray_Type, &namen_tmp, &PyArray_Type, &oldnamen_tmp, &PyArray_Type, &volume_tmp, &PyArray_Type, &centergn_tmp, &PyArray_Type, &halocentergn_tmp, &PyArray_Type, &periodicnid_tmp, &PyArray_Type, &cellnid_tmp, &PyArray_Type, &centerc_tmp, &PyArray_Type, &centerh_tmp, &PyArray_Type, &halonid_tmp, &PyArray_Type, &airDiamond_tmp, &PyArray_Type, &lambda_x_tmp, &PyArray_Type, &lambda_y_tmp, &PyArray_Type, &lambda_z_tmp, &PyArray_Type, &number_tmp, &PyArray_Type, &R_x_tmp, &PyArray_Type, &R_y_tmp, &PyArray_Type, &R_z_tmp, &PyArray_Type, &param1_tmp, &PyArray_Type, &param2_tmp, &PyArray_Type, &param3_tmp, &PyArray_Type, &shift_tmp, &PyArray_Type, &loctoglob_tmp, &PyArray_Type, &BCdirichlet_tmp, &PyArray_Type, &a_loc_tmp, &PyArray_Type, &irn_loc_tmp, &PyArray_Type, &jcn_loc_tmp, &PyArray_Type, &matrixinnerfaces_tmp, &PyArray_Type, &halofaces_tmp, &PyArray_Type, &dirichletfaces_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(cellfid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellfid = pyarray_to_ndarray(cellfid_tmp);
    }
    if (!pyarray_check(nodeidf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nodeidf = pyarray_to_ndarray(nodeidf_tmp);
    }
    if (!pyarray_check(vertexn_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        vertexn = pyarray_to_ndarray(vertexn_tmp);
    }
    if (!pyarray_check(halofid_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofid = pyarray_to_ndarray(halofid_tmp);
    }
    if (!pyarray_check(haloext_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        haloext = pyarray_to_ndarray(haloext_tmp);
    }
    if (!pyarray_check(namen_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        namen = pyarray_to_ndarray(namen_tmp);
    }
    if (!pyarray_check(oldnamen_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        oldnamen = pyarray_to_ndarray(oldnamen_tmp);
    }
    if (!pyarray_check(volume_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        volume = pyarray_to_ndarray(volume_tmp);
    }
    if (!pyarray_check(centergn_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centergn = pyarray_to_ndarray(centergn_tmp);
    }
    if (!pyarray_check(halocentergn_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        halocentergn = pyarray_to_ndarray(halocentergn_tmp);
    }
    if (!pyarray_check(periodicnid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        periodicnid = pyarray_to_ndarray(periodicnid_tmp);
    }
    if (!pyarray_check(cellnid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellnid = pyarray_to_ndarray(cellnid_tmp);
    }
    if (!pyarray_check(centerc_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerc = pyarray_to_ndarray(centerc_tmp);
    }
    if (!pyarray_check(centerh_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerh = pyarray_to_ndarray(centerh_tmp);
    }
    if (!pyarray_check(halonid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        halonid = pyarray_to_ndarray(halonid_tmp);
    }
    if (!pyarray_check(airDiamond_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        airDiamond = pyarray_to_ndarray(airDiamond_tmp);
    }
    if (!pyarray_check(lambda_x_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        lambda_x = pyarray_to_ndarray(lambda_x_tmp);
    }
    if (!pyarray_check(lambda_y_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        lambda_y = pyarray_to_ndarray(lambda_y_tmp);
    }
    if (!pyarray_check(lambda_z_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        lambda_z = pyarray_to_ndarray(lambda_z_tmp);
    }
    if (!pyarray_check(number_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        number = pyarray_to_ndarray(number_tmp);
    }
    if (!pyarray_check(R_x_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        R_x = pyarray_to_ndarray(R_x_tmp);
    }
    if (!pyarray_check(R_y_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        R_y = pyarray_to_ndarray(R_y_tmp);
    }
    if (!pyarray_check(R_z_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        R_z = pyarray_to_ndarray(R_z_tmp);
    }
    if (!pyarray_check(param1_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        param1 = pyarray_to_ndarray(param1_tmp);
    }
    if (!pyarray_check(param2_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        param2 = pyarray_to_ndarray(param2_tmp);
    }
    if (!pyarray_check(param3_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        param3 = pyarray_to_ndarray(param3_tmp);
    }
    if (!pyarray_check(shift_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        shift = pyarray_to_ndarray(shift_tmp);
    }
    if (!pyarray_check(loctoglob_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        loctoglob = pyarray_to_ndarray(loctoglob_tmp);
    }
    if (!pyarray_check(BCdirichlet_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        BCdirichlet = pyarray_to_ndarray(BCdirichlet_tmp);
    }
    if (!pyarray_check(a_loc_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        a_loc = pyarray_to_ndarray(a_loc_tmp);
    }
    if (!pyarray_check(irn_loc_tmp, NPY_INT32, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        irn_loc = pyarray_to_ndarray(irn_loc_tmp);
    }
    if (!pyarray_check(jcn_loc_tmp, NPY_INT32, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        jcn_loc = pyarray_to_ndarray(jcn_loc_tmp);
    }
    if (!pyarray_check(matrixinnerfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        matrixinnerfaces = pyarray_to_ndarray(matrixinnerfaces_tmp);
    }
    if (!pyarray_check(halofaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofaces = pyarray_to_ndarray(halofaces_tmp);
    }
    if (!pyarray_check(dirichletfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        dirichletfaces = pyarray_to_ndarray(dirichletfaces_tmp);
    }
    bind_c_get_triplet_3d(nd_ndim(&cellfid, 0), nd_ndim(&cellfid, 1), nd_data(&cellfid), nd_ndim(&nodeidf, 0), nd_ndim(&nodeidf, 1), nd_data(&nodeidf), nd_ndim(&vertexn, 0), nd_ndim(&vertexn, 1), nd_data(&vertexn), nd_ndim(&halofid, 0), nd_data(&halofid), nd_ndim(&haloext, 0), nd_ndim(&haloext, 1), nd_data(&haloext), nd_ndim(&namen, 0), nd_data(&namen), nd_ndim(&oldnamen, 0), nd_data(&oldnamen), nd_ndim(&volume, 0), nd_data(&volume), nd_ndim(&centergn, 0), nd_ndim(&centergn, 1), nd_ndim(&centergn, 2), nd_data(&centergn), nd_ndim(&halocentergn, 0), nd_ndim(&halocentergn, 1), nd_ndim(&halocentergn, 2), nd_data(&halocentergn), nd_ndim(&periodicnid, 0), nd_ndim(&periodicnid, 1), nd_data(&periodicnid), nd_ndim(&cellnid, 0), nd_ndim(&cellnid, 1), nd_data(&cellnid), nd_ndim(&centerc, 0), nd_ndim(&centerc, 1), nd_data(&centerc), nd_ndim(&centerh, 0), nd_ndim(&centerh, 1), nd_data(&centerh), nd_ndim(&halonid, 0), nd_ndim(&halonid, 1), nd_data(&halonid), nd_ndim(&airDiamond, 0), nd_data(&airDiamond), nd_ndim(&lambda_x, 0), nd_data(&lambda_x), nd_ndim(&lambda_y, 0), nd_data(&lambda_y), nd_ndim(&lambda_z, 0), nd_data(&lambda_z), nd_ndim(&number, 0), nd_data(&number), nd_ndim(&R_x, 0), nd_data(&R_x), nd_ndim(&R_y, 0), nd_data(&R_y), nd_ndim(&R_z, 0), nd_data(&R_z), nd_ndim(&param1, 0), nd_data(&param1), nd_ndim(&param2, 0), nd_data(&param2), nd_ndim(&param3, 0), nd_data(&param3), nd_ndim(&shift, 0), nd_ndim(&shift, 1), nd_data(&shift), nd_ndim(&loctoglob, 0), nd_data(&loctoglob), nd_ndim(&BCdirichlet, 0), nd_data(&BCdirichlet), nd_ndim(&a_loc, 0), nd_data(&a_loc), nd_ndim(&irn_loc, 0), nd_data(&irn_loc), nd_ndim(&jcn_loc, 0), nd_data(&jcn_loc), nd_ndim(&matrixinnerfaces, 0), nd_data(&matrixinnerfaces), nd_ndim(&halofaces, 0), nd_data(&halofaces), nd_ndim(&dirichletfaces, 0), nd_data(&dirichletfaces));
    result = Py_BuildValue("");
    free_pointer(cellfid);
    free_pointer(nodeidf);
    free_pointer(vertexn);
    free_pointer(halofid);
    free_pointer(haloext);
    free_pointer(namen);
    free_pointer(oldnamen);
    free_pointer(volume);
    free_pointer(centergn);
    free_pointer(halocentergn);
    free_pointer(periodicnid);
    free_pointer(cellnid);
    free_pointer(centerc);
    free_pointer(centerh);
    free_pointer(halonid);
    free_pointer(airDiamond);
    free_pointer(lambda_x);
    free_pointer(lambda_y);
    free_pointer(lambda_z);
    free_pointer(number);
    free_pointer(R_x);
    free_pointer(R_y);
    free_pointer(R_z);
    free_pointer(param1);
    free_pointer(param2);
    free_pointer(param3);
    free_pointer(shift);
    free_pointer(loctoglob);
    free_pointer(BCdirichlet);
    free_pointer(a_loc);
    free_pointer(irn_loc);
    free_pointer(jcn_loc);
    free_pointer(matrixinnerfaces);
    free_pointer(halofaces);
    free_pointer(dirichletfaces);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *get_rhs_loc_2d_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray cellfid = {.shape = NULL};
    t_ndarray nodeidf = {.shape = NULL};
    t_ndarray oldname = {.shape = NULL};
    t_ndarray volume = {.shape = NULL};
    t_ndarray centergn = {.shape = NULL};
    t_ndarray param1 = {.shape = NULL};
    t_ndarray param2 = {.shape = NULL};
    t_ndarray param3 = {.shape = NULL};
    t_ndarray param4 = {.shape = NULL};
    t_ndarray Pbordnode = {.shape = NULL};
    t_ndarray Pbordface = {.shape = NULL};
    t_ndarray rhs_loc = {.shape = NULL};
    t_ndarray BCdirichlet = {.shape = NULL};
    t_ndarray centergf = {.shape = NULL};
    t_ndarray matrixinnerfaces = {.shape = NULL};
    t_ndarray halofaces = {.shape = NULL};
    t_ndarray dirichletfaces = {.shape = NULL};
    PyArrayObject *cellfid_tmp;
    PyArrayObject *nodeidf_tmp;
    PyArrayObject *oldname_tmp;
    PyArrayObject *volume_tmp;
    PyArrayObject *centergn_tmp;
    PyArrayObject *param1_tmp;
    PyArrayObject *param2_tmp;
    PyArrayObject *param3_tmp;
    PyArrayObject *param4_tmp;
    PyArrayObject *Pbordnode_tmp;
    PyArrayObject *Pbordface_tmp;
    PyArrayObject *rhs_loc_tmp;
    PyArrayObject *BCdirichlet_tmp;
    PyArrayObject *centergf_tmp;
    PyArrayObject *matrixinnerfaces_tmp;
    PyArrayObject *halofaces_tmp;
    PyArrayObject *dirichletfaces_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "cellfid",
        "nodeidf",
        "oldname",
        "volume",
        "centergn",
        "param1",
        "param2",
        "param3",
        "param4",
        "Pbordnode",
        "Pbordface",
        "rhs_loc",
        "BCdirichlet",
        "centergf",
        "matrixinnerfaces",
        "halofaces",
        "dirichletfaces",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!", kwlist, &PyArray_Type, &cellfid_tmp, &PyArray_Type, &nodeidf_tmp, &PyArray_Type, &oldname_tmp, &PyArray_Type, &volume_tmp, &PyArray_Type, &centergn_tmp, &PyArray_Type, &param1_tmp, &PyArray_Type, &param2_tmp, &PyArray_Type, &param3_tmp, &PyArray_Type, &param4_tmp, &PyArray_Type, &Pbordnode_tmp, &PyArray_Type, &Pbordface_tmp, &PyArray_Type, &rhs_loc_tmp, &PyArray_Type, &BCdirichlet_tmp, &PyArray_Type, &centergf_tmp, &PyArray_Type, &matrixinnerfaces_tmp, &PyArray_Type, &halofaces_tmp, &PyArray_Type, &dirichletfaces_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(cellfid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellfid = pyarray_to_ndarray(cellfid_tmp);
    }
    if (!pyarray_check(nodeidf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nodeidf = pyarray_to_ndarray(nodeidf_tmp);
    }
    if (!pyarray_check(oldname_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        oldname = pyarray_to_ndarray(oldname_tmp);
    }
    if (!pyarray_check(volume_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        volume = pyarray_to_ndarray(volume_tmp);
    }
    if (!pyarray_check(centergn_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centergn = pyarray_to_ndarray(centergn_tmp);
    }
    if (!pyarray_check(param1_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        param1 = pyarray_to_ndarray(param1_tmp);
    }
    if (!pyarray_check(param2_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        param2 = pyarray_to_ndarray(param2_tmp);
    }
    if (!pyarray_check(param3_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        param3 = pyarray_to_ndarray(param3_tmp);
    }
    if (!pyarray_check(param4_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        param4 = pyarray_to_ndarray(param4_tmp);
    }
    if (!pyarray_check(Pbordnode_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        Pbordnode = pyarray_to_ndarray(Pbordnode_tmp);
    }
    if (!pyarray_check(Pbordface_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        Pbordface = pyarray_to_ndarray(Pbordface_tmp);
    }
    if (!pyarray_check(rhs_loc_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        rhs_loc = pyarray_to_ndarray(rhs_loc_tmp);
    }
    if (!pyarray_check(BCdirichlet_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        BCdirichlet = pyarray_to_ndarray(BCdirichlet_tmp);
    }
    if (!pyarray_check(centergf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centergf = pyarray_to_ndarray(centergf_tmp);
    }
    if (!pyarray_check(matrixinnerfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        matrixinnerfaces = pyarray_to_ndarray(matrixinnerfaces_tmp);
    }
    if (!pyarray_check(halofaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofaces = pyarray_to_ndarray(halofaces_tmp);
    }
    if (!pyarray_check(dirichletfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        dirichletfaces = pyarray_to_ndarray(dirichletfaces_tmp);
    }
    bind_c_get_rhs_loc_2d(nd_ndim(&cellfid, 0), nd_ndim(&cellfid, 1), nd_data(&cellfid), nd_ndim(&nodeidf, 0), nd_ndim(&nodeidf, 1), nd_data(&nodeidf), nd_ndim(&oldname, 0), nd_data(&oldname), nd_ndim(&volume, 0), nd_data(&volume), nd_ndim(&centergn, 0), nd_ndim(&centergn, 1), nd_ndim(&centergn, 2), nd_data(&centergn), nd_ndim(&param1, 0), nd_data(&param1), nd_ndim(&param2, 0), nd_data(&param2), nd_ndim(&param3, 0), nd_data(&param3), nd_ndim(&param4, 0), nd_data(&param4), nd_ndim(&Pbordnode, 0), nd_data(&Pbordnode), nd_ndim(&Pbordface, 0), nd_data(&Pbordface), nd_ndim(&rhs_loc, 0), nd_data(&rhs_loc), nd_ndim(&BCdirichlet, 0), nd_data(&BCdirichlet), nd_ndim(&centergf, 0), nd_ndim(&centergf, 1), nd_data(&centergf), nd_ndim(&matrixinnerfaces, 0), nd_data(&matrixinnerfaces), nd_ndim(&halofaces, 0), nd_data(&halofaces), nd_ndim(&dirichletfaces, 0), nd_data(&dirichletfaces));
    result = Py_BuildValue("");
    free_pointer(cellfid);
    free_pointer(nodeidf);
    free_pointer(oldname);
    free_pointer(volume);
    free_pointer(centergn);
    free_pointer(param1);
    free_pointer(param2);
    free_pointer(param3);
    free_pointer(param4);
    free_pointer(Pbordnode);
    free_pointer(Pbordface);
    free_pointer(rhs_loc);
    free_pointer(BCdirichlet);
    free_pointer(centergf);
    free_pointer(matrixinnerfaces);
    free_pointer(halofaces);
    free_pointer(dirichletfaces);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *get_rhs_glob_2d_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray cellfid = {.shape = NULL};
    t_ndarray nodeidf = {.shape = NULL};
    t_ndarray oldname = {.shape = NULL};
    t_ndarray volume = {.shape = NULL};
    t_ndarray centergn = {.shape = NULL};
    t_ndarray loctoglob = {.shape = NULL};
    t_ndarray param1 = {.shape = NULL};
    t_ndarray param2 = {.shape = NULL};
    t_ndarray param3 = {.shape = NULL};
    t_ndarray param4 = {.shape = NULL};
    t_ndarray Pbordnode = {.shape = NULL};
    t_ndarray Pbordface = {.shape = NULL};
    t_ndarray rhs = {.shape = NULL};
    t_ndarray BCdirichlet = {.shape = NULL};
    t_ndarray centergf = {.shape = NULL};
    t_ndarray matrixinnerfaces = {.shape = NULL};
    t_ndarray halofaces = {.shape = NULL};
    t_ndarray dirichletfaces = {.shape = NULL};
    PyArrayObject *cellfid_tmp;
    PyArrayObject *nodeidf_tmp;
    PyArrayObject *oldname_tmp;
    PyArrayObject *volume_tmp;
    PyArrayObject *centergn_tmp;
    PyArrayObject *loctoglob_tmp;
    PyArrayObject *param1_tmp;
    PyArrayObject *param2_tmp;
    PyArrayObject *param3_tmp;
    PyArrayObject *param4_tmp;
    PyArrayObject *Pbordnode_tmp;
    PyArrayObject *Pbordface_tmp;
    PyArrayObject *rhs_tmp;
    PyArrayObject *BCdirichlet_tmp;
    PyArrayObject *centergf_tmp;
    PyArrayObject *matrixinnerfaces_tmp;
    PyArrayObject *halofaces_tmp;
    PyArrayObject *dirichletfaces_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "cellfid",
        "nodeidf",
        "oldname",
        "volume",
        "centergn",
        "loctoglob",
        "param1",
        "param2",
        "param3",
        "param4",
        "Pbordnode",
        "Pbordface",
        "rhs",
        "BCdirichlet",
        "centergf",
        "matrixinnerfaces",
        "halofaces",
        "dirichletfaces",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!", kwlist, &PyArray_Type, &cellfid_tmp, &PyArray_Type, &nodeidf_tmp, &PyArray_Type, &oldname_tmp, &PyArray_Type, &volume_tmp, &PyArray_Type, &centergn_tmp, &PyArray_Type, &loctoglob_tmp, &PyArray_Type, &param1_tmp, &PyArray_Type, &param2_tmp, &PyArray_Type, &param3_tmp, &PyArray_Type, &param4_tmp, &PyArray_Type, &Pbordnode_tmp, &PyArray_Type, &Pbordface_tmp, &PyArray_Type, &rhs_tmp, &PyArray_Type, &BCdirichlet_tmp, &PyArray_Type, &centergf_tmp, &PyArray_Type, &matrixinnerfaces_tmp, &PyArray_Type, &halofaces_tmp, &PyArray_Type, &dirichletfaces_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(cellfid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellfid = pyarray_to_ndarray(cellfid_tmp);
    }
    if (!pyarray_check(nodeidf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nodeidf = pyarray_to_ndarray(nodeidf_tmp);
    }
    if (!pyarray_check(oldname_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        oldname = pyarray_to_ndarray(oldname_tmp);
    }
    if (!pyarray_check(volume_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        volume = pyarray_to_ndarray(volume_tmp);
    }
    if (!pyarray_check(centergn_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centergn = pyarray_to_ndarray(centergn_tmp);
    }
    if (!pyarray_check(loctoglob_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        loctoglob = pyarray_to_ndarray(loctoglob_tmp);
    }
    if (!pyarray_check(param1_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        param1 = pyarray_to_ndarray(param1_tmp);
    }
    if (!pyarray_check(param2_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        param2 = pyarray_to_ndarray(param2_tmp);
    }
    if (!pyarray_check(param3_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        param3 = pyarray_to_ndarray(param3_tmp);
    }
    if (!pyarray_check(param4_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        param4 = pyarray_to_ndarray(param4_tmp);
    }
    if (!pyarray_check(Pbordnode_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        Pbordnode = pyarray_to_ndarray(Pbordnode_tmp);
    }
    if (!pyarray_check(Pbordface_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        Pbordface = pyarray_to_ndarray(Pbordface_tmp);
    }
    if (!pyarray_check(rhs_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        rhs = pyarray_to_ndarray(rhs_tmp);
    }
    if (!pyarray_check(BCdirichlet_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        BCdirichlet = pyarray_to_ndarray(BCdirichlet_tmp);
    }
    if (!pyarray_check(centergf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centergf = pyarray_to_ndarray(centergf_tmp);
    }
    if (!pyarray_check(matrixinnerfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        matrixinnerfaces = pyarray_to_ndarray(matrixinnerfaces_tmp);
    }
    if (!pyarray_check(halofaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofaces = pyarray_to_ndarray(halofaces_tmp);
    }
    if (!pyarray_check(dirichletfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        dirichletfaces = pyarray_to_ndarray(dirichletfaces_tmp);
    }
    bind_c_get_rhs_glob_2d(nd_ndim(&cellfid, 0), nd_ndim(&cellfid, 1), nd_data(&cellfid), nd_ndim(&nodeidf, 0), nd_ndim(&nodeidf, 1), nd_data(&nodeidf), nd_ndim(&oldname, 0), nd_data(&oldname), nd_ndim(&volume, 0), nd_data(&volume), nd_ndim(&centergn, 0), nd_ndim(&centergn, 1), nd_ndim(&centergn, 2), nd_data(&centergn), nd_ndim(&loctoglob, 0), nd_data(&loctoglob), nd_ndim(&param1, 0), nd_data(&param1), nd_ndim(&param2, 0), nd_data(&param2), nd_ndim(&param3, 0), nd_data(&param3), nd_ndim(&param4, 0), nd_data(&param4), nd_ndim(&Pbordnode, 0), nd_data(&Pbordnode), nd_ndim(&Pbordface, 0), nd_data(&Pbordface), nd_ndim(&rhs, 0), nd_data(&rhs), nd_ndim(&BCdirichlet, 0), nd_data(&BCdirichlet), nd_ndim(&centergf, 0), nd_ndim(&centergf, 1), nd_data(&centergf), nd_ndim(&matrixinnerfaces, 0), nd_data(&matrixinnerfaces), nd_ndim(&halofaces, 0), nd_data(&halofaces), nd_ndim(&dirichletfaces, 0), nd_data(&dirichletfaces));
    result = Py_BuildValue("");
    free_pointer(cellfid);
    free_pointer(nodeidf);
    free_pointer(oldname);
    free_pointer(volume);
    free_pointer(centergn);
    free_pointer(loctoglob);
    free_pointer(param1);
    free_pointer(param2);
    free_pointer(param3);
    free_pointer(param4);
    free_pointer(Pbordnode);
    free_pointer(Pbordface);
    free_pointer(rhs);
    free_pointer(BCdirichlet);
    free_pointer(centergf);
    free_pointer(matrixinnerfaces);
    free_pointer(halofaces);
    free_pointer(dirichletfaces);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *get_rhs_loc_3d_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray cellfid = {.shape = NULL};
    t_ndarray nodeidf = {.shape = NULL};
    t_ndarray oldname = {.shape = NULL};
    t_ndarray volume = {.shape = NULL};
    t_ndarray centergn = {.shape = NULL};
    t_ndarray param1 = {.shape = NULL};
    t_ndarray param2 = {.shape = NULL};
    t_ndarray param3 = {.shape = NULL};
    t_ndarray Pbordnode = {.shape = NULL};
    t_ndarray Pbordface = {.shape = NULL};
    t_ndarray rhs_loc = {.shape = NULL};
    t_ndarray BCdirichlet = {.shape = NULL};
    t_ndarray matrixinnerfaces = {.shape = NULL};
    t_ndarray halofaces = {.shape = NULL};
    t_ndarray dirichletfaces = {.shape = NULL};
    PyArrayObject *cellfid_tmp;
    PyArrayObject *nodeidf_tmp;
    PyArrayObject *oldname_tmp;
    PyArrayObject *volume_tmp;
    PyArrayObject *centergn_tmp;
    PyArrayObject *param1_tmp;
    PyArrayObject *param2_tmp;
    PyArrayObject *param3_tmp;
    PyArrayObject *Pbordnode_tmp;
    PyArrayObject *Pbordface_tmp;
    PyArrayObject *rhs_loc_tmp;
    PyArrayObject *BCdirichlet_tmp;
    PyArrayObject *matrixinnerfaces_tmp;
    PyArrayObject *halofaces_tmp;
    PyArrayObject *dirichletfaces_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "cellfid",
        "nodeidf",
        "oldname",
        "volume",
        "centergn",
        "param1",
        "param2",
        "param3",
        "Pbordnode",
        "Pbordface",
        "rhs_loc",
        "BCdirichlet",
        "matrixinnerfaces",
        "halofaces",
        "dirichletfaces",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!", kwlist, &PyArray_Type, &cellfid_tmp, &PyArray_Type, &nodeidf_tmp, &PyArray_Type, &oldname_tmp, &PyArray_Type, &volume_tmp, &PyArray_Type, &centergn_tmp, &PyArray_Type, &param1_tmp, &PyArray_Type, &param2_tmp, &PyArray_Type, &param3_tmp, &PyArray_Type, &Pbordnode_tmp, &PyArray_Type, &Pbordface_tmp, &PyArray_Type, &rhs_loc_tmp, &PyArray_Type, &BCdirichlet_tmp, &PyArray_Type, &matrixinnerfaces_tmp, &PyArray_Type, &halofaces_tmp, &PyArray_Type, &dirichletfaces_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(cellfid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellfid = pyarray_to_ndarray(cellfid_tmp);
    }
    if (!pyarray_check(nodeidf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nodeidf = pyarray_to_ndarray(nodeidf_tmp);
    }
    if (!pyarray_check(oldname_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        oldname = pyarray_to_ndarray(oldname_tmp);
    }
    if (!pyarray_check(volume_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        volume = pyarray_to_ndarray(volume_tmp);
    }
    if (!pyarray_check(centergn_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centergn = pyarray_to_ndarray(centergn_tmp);
    }
    if (!pyarray_check(param1_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        param1 = pyarray_to_ndarray(param1_tmp);
    }
    if (!pyarray_check(param2_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        param2 = pyarray_to_ndarray(param2_tmp);
    }
    if (!pyarray_check(param3_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        param3 = pyarray_to_ndarray(param3_tmp);
    }
    if (!pyarray_check(Pbordnode_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        Pbordnode = pyarray_to_ndarray(Pbordnode_tmp);
    }
    if (!pyarray_check(Pbordface_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        Pbordface = pyarray_to_ndarray(Pbordface_tmp);
    }
    if (!pyarray_check(rhs_loc_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        rhs_loc = pyarray_to_ndarray(rhs_loc_tmp);
    }
    if (!pyarray_check(BCdirichlet_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        BCdirichlet = pyarray_to_ndarray(BCdirichlet_tmp);
    }
    if (!pyarray_check(matrixinnerfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        matrixinnerfaces = pyarray_to_ndarray(matrixinnerfaces_tmp);
    }
    if (!pyarray_check(halofaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofaces = pyarray_to_ndarray(halofaces_tmp);
    }
    if (!pyarray_check(dirichletfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        dirichletfaces = pyarray_to_ndarray(dirichletfaces_tmp);
    }
    bind_c_get_rhs_loc_3d(nd_ndim(&cellfid, 0), nd_ndim(&cellfid, 1), nd_data(&cellfid), nd_ndim(&nodeidf, 0), nd_ndim(&nodeidf, 1), nd_data(&nodeidf), nd_ndim(&oldname, 0), nd_data(&oldname), nd_ndim(&volume, 0), nd_data(&volume), nd_ndim(&centergn, 0), nd_ndim(&centergn, 1), nd_ndim(&centergn, 2), nd_data(&centergn), nd_ndim(&param1, 0), nd_data(&param1), nd_ndim(&param2, 0), nd_data(&param2), nd_ndim(&param3, 0), nd_data(&param3), nd_ndim(&Pbordnode, 0), nd_data(&Pbordnode), nd_ndim(&Pbordface, 0), nd_data(&Pbordface), nd_ndim(&rhs_loc, 0), nd_data(&rhs_loc), nd_ndim(&BCdirichlet, 0), nd_data(&BCdirichlet), nd_ndim(&matrixinnerfaces, 0), nd_data(&matrixinnerfaces), nd_ndim(&halofaces, 0), nd_data(&halofaces), nd_ndim(&dirichletfaces, 0), nd_data(&dirichletfaces));
    result = Py_BuildValue("");
    free_pointer(cellfid);
    free_pointer(nodeidf);
    free_pointer(oldname);
    free_pointer(volume);
    free_pointer(centergn);
    free_pointer(param1);
    free_pointer(param2);
    free_pointer(param3);
    free_pointer(Pbordnode);
    free_pointer(Pbordface);
    free_pointer(rhs_loc);
    free_pointer(BCdirichlet);
    free_pointer(matrixinnerfaces);
    free_pointer(halofaces);
    free_pointer(dirichletfaces);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *get_rhs_glob_3d_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray cellfid = {.shape = NULL};
    t_ndarray nodeidf = {.shape = NULL};
    t_ndarray oldname = {.shape = NULL};
    t_ndarray volume = {.shape = NULL};
    t_ndarray centergn = {.shape = NULL};
    t_ndarray loctoglob = {.shape = NULL};
    t_ndarray param1 = {.shape = NULL};
    t_ndarray param2 = {.shape = NULL};
    t_ndarray param3 = {.shape = NULL};
    t_ndarray Pbordnode = {.shape = NULL};
    t_ndarray Pbordface = {.shape = NULL};
    t_ndarray rhs = {.shape = NULL};
    t_ndarray BCdirichlet = {.shape = NULL};
    t_ndarray matrixinnerfaces = {.shape = NULL};
    t_ndarray halofaces = {.shape = NULL};
    t_ndarray dirichletfaces = {.shape = NULL};
    PyArrayObject *cellfid_tmp;
    PyArrayObject *nodeidf_tmp;
    PyArrayObject *oldname_tmp;
    PyArrayObject *volume_tmp;
    PyArrayObject *centergn_tmp;
    PyArrayObject *loctoglob_tmp;
    PyArrayObject *param1_tmp;
    PyArrayObject *param2_tmp;
    PyArrayObject *param3_tmp;
    PyArrayObject *Pbordnode_tmp;
    PyArrayObject *Pbordface_tmp;
    PyArrayObject *rhs_tmp;
    PyArrayObject *BCdirichlet_tmp;
    PyArrayObject *matrixinnerfaces_tmp;
    PyArrayObject *halofaces_tmp;
    PyArrayObject *dirichletfaces_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "cellfid",
        "nodeidf",
        "oldname",
        "volume",
        "centergn",
        "loctoglob",
        "param1",
        "param2",
        "param3",
        "Pbordnode",
        "Pbordface",
        "rhs",
        "BCdirichlet",
        "matrixinnerfaces",
        "halofaces",
        "dirichletfaces",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!", kwlist, &PyArray_Type, &cellfid_tmp, &PyArray_Type, &nodeidf_tmp, &PyArray_Type, &oldname_tmp, &PyArray_Type, &volume_tmp, &PyArray_Type, &centergn_tmp, &PyArray_Type, &loctoglob_tmp, &PyArray_Type, &param1_tmp, &PyArray_Type, &param2_tmp, &PyArray_Type, &param3_tmp, &PyArray_Type, &Pbordnode_tmp, &PyArray_Type, &Pbordface_tmp, &PyArray_Type, &rhs_tmp, &PyArray_Type, &BCdirichlet_tmp, &PyArray_Type, &matrixinnerfaces_tmp, &PyArray_Type, &halofaces_tmp, &PyArray_Type, &dirichletfaces_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(cellfid_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellfid = pyarray_to_ndarray(cellfid_tmp);
    }
    if (!pyarray_check(nodeidf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nodeidf = pyarray_to_ndarray(nodeidf_tmp);
    }
    if (!pyarray_check(oldname_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        oldname = pyarray_to_ndarray(oldname_tmp);
    }
    if (!pyarray_check(volume_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        volume = pyarray_to_ndarray(volume_tmp);
    }
    if (!pyarray_check(centergn_tmp, NPY_DOUBLE, 3, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centergn = pyarray_to_ndarray(centergn_tmp);
    }
    if (!pyarray_check(loctoglob_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        loctoglob = pyarray_to_ndarray(loctoglob_tmp);
    }
    if (!pyarray_check(param1_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        param1 = pyarray_to_ndarray(param1_tmp);
    }
    if (!pyarray_check(param2_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        param2 = pyarray_to_ndarray(param2_tmp);
    }
    if (!pyarray_check(param3_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        param3 = pyarray_to_ndarray(param3_tmp);
    }
    if (!pyarray_check(Pbordnode_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        Pbordnode = pyarray_to_ndarray(Pbordnode_tmp);
    }
    if (!pyarray_check(Pbordface_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        Pbordface = pyarray_to_ndarray(Pbordface_tmp);
    }
    if (!pyarray_check(rhs_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        rhs = pyarray_to_ndarray(rhs_tmp);
    }
    if (!pyarray_check(BCdirichlet_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        BCdirichlet = pyarray_to_ndarray(BCdirichlet_tmp);
    }
    if (!pyarray_check(matrixinnerfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        matrixinnerfaces = pyarray_to_ndarray(matrixinnerfaces_tmp);
    }
    if (!pyarray_check(halofaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofaces = pyarray_to_ndarray(halofaces_tmp);
    }
    if (!pyarray_check(dirichletfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        dirichletfaces = pyarray_to_ndarray(dirichletfaces_tmp);
    }
    bind_c_get_rhs_glob_3d(nd_ndim(&cellfid, 0), nd_ndim(&cellfid, 1), nd_data(&cellfid), nd_ndim(&nodeidf, 0), nd_ndim(&nodeidf, 1), nd_data(&nodeidf), nd_ndim(&oldname, 0), nd_data(&oldname), nd_ndim(&volume, 0), nd_data(&volume), nd_ndim(&centergn, 0), nd_ndim(&centergn, 1), nd_ndim(&centergn, 2), nd_data(&centergn), nd_ndim(&loctoglob, 0), nd_data(&loctoglob), nd_ndim(&param1, 0), nd_data(&param1), nd_ndim(&param2, 0), nd_data(&param2), nd_ndim(&param3, 0), nd_data(&param3), nd_ndim(&Pbordnode, 0), nd_data(&Pbordnode), nd_ndim(&Pbordface, 0), nd_data(&Pbordface), nd_ndim(&rhs, 0), nd_data(&rhs), nd_ndim(&BCdirichlet, 0), nd_data(&BCdirichlet), nd_ndim(&matrixinnerfaces, 0), nd_data(&matrixinnerfaces), nd_ndim(&halofaces, 0), nd_data(&halofaces), nd_ndim(&dirichletfaces, 0), nd_data(&dirichletfaces));
    result = Py_BuildValue("");
    free_pointer(cellfid);
    free_pointer(nodeidf);
    free_pointer(oldname);
    free_pointer(volume);
    free_pointer(centergn);
    free_pointer(loctoglob);
    free_pointer(param1);
    free_pointer(param2);
    free_pointer(param3);
    free_pointer(Pbordnode);
    free_pointer(Pbordface);
    free_pointer(rhs);
    free_pointer(BCdirichlet);
    free_pointer(matrixinnerfaces);
    free_pointer(halofaces);
    free_pointer(dirichletfaces);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *compute_P_gradient_2d_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray P_c = {.shape = NULL};
    t_ndarray P_ghost = {.shape = NULL};
    t_ndarray P_halo = {.shape = NULL};
    t_ndarray P_node = {.shape = NULL};
    t_ndarray cellidf = {.shape = NULL};
    t_ndarray nodeidf = {.shape = NULL};
    t_ndarray centergf = {.shape = NULL};
    t_ndarray namef = {.shape = NULL};
    t_ndarray halofid = {.shape = NULL};
    t_ndarray centerc = {.shape = NULL};
    t_ndarray centerh = {.shape = NULL};
    t_ndarray oldname = {.shape = NULL};
    t_ndarray airDiamond = {.shape = NULL};
    t_ndarray f_1 = {.shape = NULL};
    t_ndarray f_2 = {.shape = NULL};
    t_ndarray f_3 = {.shape = NULL};
    t_ndarray f_4 = {.shape = NULL};
    t_ndarray normalf = {.shape = NULL};
    t_ndarray shift = {.shape = NULL};
    t_ndarray Pbordnode = {.shape = NULL};
    t_ndarray Pbordface = {.shape = NULL};
    t_ndarray Px_face = {.shape = NULL};
    t_ndarray Py_face = {.shape = NULL};
    t_ndarray Pz_face = {.shape = NULL};
    t_ndarray BCdirichlet = {.shape = NULL};
    t_ndarray innerfaces = {.shape = NULL};
    t_ndarray halofaces = {.shape = NULL};
    t_ndarray neumannfaces = {.shape = NULL};
    t_ndarray dirichletfaces = {.shape = NULL};
    t_ndarray periodicfaces = {.shape = NULL};
    PyArrayObject *P_c_tmp;
    PyArrayObject *P_ghost_tmp;
    PyArrayObject *P_halo_tmp;
    PyArrayObject *P_node_tmp;
    PyArrayObject *cellidf_tmp;
    PyArrayObject *nodeidf_tmp;
    PyArrayObject *centergf_tmp;
    PyArrayObject *namef_tmp;
    PyArrayObject *halofid_tmp;
    PyArrayObject *centerc_tmp;
    PyArrayObject *centerh_tmp;
    PyArrayObject *oldname_tmp;
    PyArrayObject *airDiamond_tmp;
    PyArrayObject *f_1_tmp;
    PyArrayObject *f_2_tmp;
    PyArrayObject *f_3_tmp;
    PyArrayObject *f_4_tmp;
    PyArrayObject *normalf_tmp;
    PyArrayObject *shift_tmp;
    PyArrayObject *Pbordnode_tmp;
    PyArrayObject *Pbordface_tmp;
    PyArrayObject *Px_face_tmp;
    PyArrayObject *Py_face_tmp;
    PyArrayObject *Pz_face_tmp;
    PyArrayObject *BCdirichlet_tmp;
    PyArrayObject *innerfaces_tmp;
    PyArrayObject *halofaces_tmp;
    PyArrayObject *neumannfaces_tmp;
    PyArrayObject *dirichletfaces_tmp;
    PyArrayObject *periodicfaces_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "P_c",
        "P_ghost",
        "P_halo",
        "P_node",
        "cellidf",
        "nodeidf",
        "centergf",
        "namef",
        "halofid",
        "centerc",
        "centerh",
        "oldname",
        "airDiamond",
        "f_1",
        "f_2",
        "f_3",
        "f_4",
        "normalf",
        "shift",
        "Pbordnode",
        "Pbordface",
        "Px_face",
        "Py_face",
        "Pz_face",
        "BCdirichlet",
        "innerfaces",
        "halofaces",
        "neumannfaces",
        "dirichletfaces",
        "periodicfaces",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!", kwlist, &PyArray_Type, &P_c_tmp, &PyArray_Type, &P_ghost_tmp, &PyArray_Type, &P_halo_tmp, &PyArray_Type, &P_node_tmp, &PyArray_Type, &cellidf_tmp, &PyArray_Type, &nodeidf_tmp, &PyArray_Type, &centergf_tmp, &PyArray_Type, &namef_tmp, &PyArray_Type, &halofid_tmp, &PyArray_Type, &centerc_tmp, &PyArray_Type, &centerh_tmp, &PyArray_Type, &oldname_tmp, &PyArray_Type, &airDiamond_tmp, &PyArray_Type, &f_1_tmp, &PyArray_Type, &f_2_tmp, &PyArray_Type, &f_3_tmp, &PyArray_Type, &f_4_tmp, &PyArray_Type, &normalf_tmp, &PyArray_Type, &shift_tmp, &PyArray_Type, &Pbordnode_tmp, &PyArray_Type, &Pbordface_tmp, &PyArray_Type, &Px_face_tmp, &PyArray_Type, &Py_face_tmp, &PyArray_Type, &Pz_face_tmp, &PyArray_Type, &BCdirichlet_tmp, &PyArray_Type, &innerfaces_tmp, &PyArray_Type, &halofaces_tmp, &PyArray_Type, &neumannfaces_tmp, &PyArray_Type, &dirichletfaces_tmp, &PyArray_Type, &periodicfaces_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(P_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        P_c = pyarray_to_ndarray(P_c_tmp);
    }
    if (!pyarray_check(P_ghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        P_ghost = pyarray_to_ndarray(P_ghost_tmp);
    }
    if (!pyarray_check(P_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        P_halo = pyarray_to_ndarray(P_halo_tmp);
    }
    if (!pyarray_check(P_node_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        P_node = pyarray_to_ndarray(P_node_tmp);
    }
    if (!pyarray_check(cellidf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellidf = pyarray_to_ndarray(cellidf_tmp);
    }
    if (!pyarray_check(nodeidf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nodeidf = pyarray_to_ndarray(nodeidf_tmp);
    }
    if (!pyarray_check(centergf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centergf = pyarray_to_ndarray(centergf_tmp);
    }
    if (!pyarray_check(namef_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        namef = pyarray_to_ndarray(namef_tmp);
    }
    if (!pyarray_check(halofid_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofid = pyarray_to_ndarray(halofid_tmp);
    }
    if (!pyarray_check(centerc_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerc = pyarray_to_ndarray(centerc_tmp);
    }
    if (!pyarray_check(centerh_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerh = pyarray_to_ndarray(centerh_tmp);
    }
    if (!pyarray_check(oldname_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        oldname = pyarray_to_ndarray(oldname_tmp);
    }
    if (!pyarray_check(airDiamond_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        airDiamond = pyarray_to_ndarray(airDiamond_tmp);
    }
    if (!pyarray_check(f_1_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        f_1 = pyarray_to_ndarray(f_1_tmp);
    }
    if (!pyarray_check(f_2_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        f_2 = pyarray_to_ndarray(f_2_tmp);
    }
    if (!pyarray_check(f_3_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        f_3 = pyarray_to_ndarray(f_3_tmp);
    }
    if (!pyarray_check(f_4_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        f_4 = pyarray_to_ndarray(f_4_tmp);
    }
    if (!pyarray_check(normalf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        normalf = pyarray_to_ndarray(normalf_tmp);
    }
    if (!pyarray_check(shift_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        shift = pyarray_to_ndarray(shift_tmp);
    }
    if (!pyarray_check(Pbordnode_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        Pbordnode = pyarray_to_ndarray(Pbordnode_tmp);
    }
    if (!pyarray_check(Pbordface_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        Pbordface = pyarray_to_ndarray(Pbordface_tmp);
    }
    if (!pyarray_check(Px_face_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        Px_face = pyarray_to_ndarray(Px_face_tmp);
    }
    if (!pyarray_check(Py_face_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        Py_face = pyarray_to_ndarray(Py_face_tmp);
    }
    if (!pyarray_check(Pz_face_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        Pz_face = pyarray_to_ndarray(Pz_face_tmp);
    }
    if (!pyarray_check(BCdirichlet_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        BCdirichlet = pyarray_to_ndarray(BCdirichlet_tmp);
    }
    if (!pyarray_check(innerfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        innerfaces = pyarray_to_ndarray(innerfaces_tmp);
    }
    if (!pyarray_check(halofaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofaces = pyarray_to_ndarray(halofaces_tmp);
    }
    if (!pyarray_check(neumannfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        neumannfaces = pyarray_to_ndarray(neumannfaces_tmp);
    }
    if (!pyarray_check(dirichletfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        dirichletfaces = pyarray_to_ndarray(dirichletfaces_tmp);
    }
    if (!pyarray_check(periodicfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        periodicfaces = pyarray_to_ndarray(periodicfaces_tmp);
    }
    bind_c_compute_p_gradient_2d(nd_ndim(&P_c, 0), nd_data(&P_c), nd_ndim(&P_ghost, 0), nd_data(&P_ghost), nd_ndim(&P_halo, 0), nd_data(&P_halo), nd_ndim(&P_node, 0), nd_data(&P_node), nd_ndim(&cellidf, 0), nd_ndim(&cellidf, 1), nd_data(&cellidf), nd_ndim(&nodeidf, 0), nd_ndim(&nodeidf, 1), nd_data(&nodeidf), nd_ndim(&centergf, 0), nd_ndim(&centergf, 1), nd_data(&centergf), nd_ndim(&namef, 0), nd_data(&namef), nd_ndim(&halofid, 0), nd_data(&halofid), nd_ndim(&centerc, 0), nd_ndim(&centerc, 1), nd_data(&centerc), nd_ndim(&centerh, 0), nd_ndim(&centerh, 1), nd_data(&centerh), nd_ndim(&oldname, 0), nd_data(&oldname), nd_ndim(&airDiamond, 0), nd_data(&airDiamond), nd_ndim(&f_1, 0), nd_ndim(&f_1, 1), nd_data(&f_1), nd_ndim(&f_2, 0), nd_ndim(&f_2, 1), nd_data(&f_2), nd_ndim(&f_3, 0), nd_ndim(&f_3, 1), nd_data(&f_3), nd_ndim(&f_4, 0), nd_ndim(&f_4, 1), nd_data(&f_4), nd_ndim(&normalf, 0), nd_ndim(&normalf, 1), nd_data(&normalf), nd_ndim(&shift, 0), nd_ndim(&shift, 1), nd_data(&shift), nd_ndim(&Pbordnode, 0), nd_data(&Pbordnode), nd_ndim(&Pbordface, 0), nd_data(&Pbordface), nd_ndim(&Px_face, 0), nd_data(&Px_face), nd_ndim(&Py_face, 0), nd_data(&Py_face), nd_ndim(&Pz_face, 0), nd_data(&Pz_face), nd_ndim(&BCdirichlet, 0), nd_data(&BCdirichlet), nd_ndim(&innerfaces, 0), nd_data(&innerfaces), nd_ndim(&halofaces, 0), nd_data(&halofaces), nd_ndim(&neumannfaces, 0), nd_data(&neumannfaces), nd_ndim(&dirichletfaces, 0), nd_data(&dirichletfaces), nd_ndim(&periodicfaces, 0), nd_data(&periodicfaces));
    result = Py_BuildValue("");
    free_pointer(P_c);
    free_pointer(P_ghost);
    free_pointer(P_halo);
    free_pointer(P_node);
    free_pointer(cellidf);
    free_pointer(nodeidf);
    free_pointer(centergf);
    free_pointer(namef);
    free_pointer(halofid);
    free_pointer(centerc);
    free_pointer(centerh);
    free_pointer(oldname);
    free_pointer(airDiamond);
    free_pointer(f_1);
    free_pointer(f_2);
    free_pointer(f_3);
    free_pointer(f_4);
    free_pointer(normalf);
    free_pointer(shift);
    free_pointer(Pbordnode);
    free_pointer(Pbordface);
    free_pointer(Px_face);
    free_pointer(Py_face);
    free_pointer(Pz_face);
    free_pointer(BCdirichlet);
    free_pointer(innerfaces);
    free_pointer(halofaces);
    free_pointer(neumannfaces);
    free_pointer(dirichletfaces);
    free_pointer(periodicfaces);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *compute_P_gradient_3d_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray val_c = {.shape = NULL};
    t_ndarray v_ghost = {.shape = NULL};
    t_ndarray v_halo = {.shape = NULL};
    t_ndarray v_node = {.shape = NULL};
    t_ndarray cellidf = {.shape = NULL};
    t_ndarray nodeidf = {.shape = NULL};
    t_ndarray centergf = {.shape = NULL};
    t_ndarray namef = {.shape = NULL};
    t_ndarray halofid = {.shape = NULL};
    t_ndarray centerc = {.shape = NULL};
    t_ndarray centerh = {.shape = NULL};
    t_ndarray oldname = {.shape = NULL};
    t_ndarray airDiamond = {.shape = NULL};
    t_ndarray n1 = {.shape = NULL};
    t_ndarray n2 = {.shape = NULL};
    t_ndarray n3 = {.shape = NULL};
    t_ndarray n4 = {.shape = NULL};
    t_ndarray normalf = {.shape = NULL};
    t_ndarray shift = {.shape = NULL};
    t_ndarray Pbordnode = {.shape = NULL};
    t_ndarray Pbordface = {.shape = NULL};
    t_ndarray Px_face = {.shape = NULL};
    t_ndarray Py_face = {.shape = NULL};
    t_ndarray Pz_face = {.shape = NULL};
    t_ndarray BCdirichlet = {.shape = NULL};
    t_ndarray innerfaces = {.shape = NULL};
    t_ndarray halofaces = {.shape = NULL};
    t_ndarray neumannfaces = {.shape = NULL};
    t_ndarray dirichletfaces = {.shape = NULL};
    t_ndarray periodicfaces = {.shape = NULL};
    PyArrayObject *val_c_tmp;
    PyArrayObject *v_ghost_tmp;
    PyArrayObject *v_halo_tmp;
    PyArrayObject *v_node_tmp;
    PyArrayObject *cellidf_tmp;
    PyArrayObject *nodeidf_tmp;
    PyArrayObject *centergf_tmp;
    PyArrayObject *namef_tmp;
    PyArrayObject *halofid_tmp;
    PyArrayObject *centerc_tmp;
    PyArrayObject *centerh_tmp;
    PyArrayObject *oldname_tmp;
    PyArrayObject *airDiamond_tmp;
    PyArrayObject *n1_tmp;
    PyArrayObject *n2_tmp;
    PyArrayObject *n3_tmp;
    PyArrayObject *n4_tmp;
    PyArrayObject *normalf_tmp;
    PyArrayObject *shift_tmp;
    PyArrayObject *Pbordnode_tmp;
    PyArrayObject *Pbordface_tmp;
    PyArrayObject *Px_face_tmp;
    PyArrayObject *Py_face_tmp;
    PyArrayObject *Pz_face_tmp;
    PyArrayObject *BCdirichlet_tmp;
    PyArrayObject *innerfaces_tmp;
    PyArrayObject *halofaces_tmp;
    PyArrayObject *neumannfaces_tmp;
    PyArrayObject *dirichletfaces_tmp;
    PyArrayObject *periodicfaces_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "val_c",
        "v_ghost",
        "v_halo",
        "v_node",
        "cellidf",
        "nodeidf",
        "centergf",
        "namef",
        "halofid",
        "centerc",
        "centerh",
        "oldname",
        "airDiamond",
        "n1",
        "n2",
        "n3",
        "n4",
        "normalf",
        "shift",
        "Pbordnode",
        "Pbordface",
        "Px_face",
        "Py_face",
        "Pz_face",
        "BCdirichlet",
        "innerfaces",
        "halofaces",
        "neumannfaces",
        "dirichletfaces",
        "periodicfaces",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!O!", kwlist, &PyArray_Type, &val_c_tmp, &PyArray_Type, &v_ghost_tmp, &PyArray_Type, &v_halo_tmp, &PyArray_Type, &v_node_tmp, &PyArray_Type, &cellidf_tmp, &PyArray_Type, &nodeidf_tmp, &PyArray_Type, &centergf_tmp, &PyArray_Type, &namef_tmp, &PyArray_Type, &halofid_tmp, &PyArray_Type, &centerc_tmp, &PyArray_Type, &centerh_tmp, &PyArray_Type, &oldname_tmp, &PyArray_Type, &airDiamond_tmp, &PyArray_Type, &n1_tmp, &PyArray_Type, &n2_tmp, &PyArray_Type, &n3_tmp, &PyArray_Type, &n4_tmp, &PyArray_Type, &normalf_tmp, &PyArray_Type, &shift_tmp, &PyArray_Type, &Pbordnode_tmp, &PyArray_Type, &Pbordface_tmp, &PyArray_Type, &Px_face_tmp, &PyArray_Type, &Py_face_tmp, &PyArray_Type, &Pz_face_tmp, &PyArray_Type, &BCdirichlet_tmp, &PyArray_Type, &innerfaces_tmp, &PyArray_Type, &halofaces_tmp, &PyArray_Type, &neumannfaces_tmp, &PyArray_Type, &dirichletfaces_tmp, &PyArray_Type, &periodicfaces_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(val_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        val_c = pyarray_to_ndarray(val_c_tmp);
    }
    if (!pyarray_check(v_ghost_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        v_ghost = pyarray_to_ndarray(v_ghost_tmp);
    }
    if (!pyarray_check(v_halo_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        v_halo = pyarray_to_ndarray(v_halo_tmp);
    }
    if (!pyarray_check(v_node_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        v_node = pyarray_to_ndarray(v_node_tmp);
    }
    if (!pyarray_check(cellidf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        cellidf = pyarray_to_ndarray(cellidf_tmp);
    }
    if (!pyarray_check(nodeidf_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        nodeidf = pyarray_to_ndarray(nodeidf_tmp);
    }
    if (!pyarray_check(centergf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centergf = pyarray_to_ndarray(centergf_tmp);
    }
    if (!pyarray_check(namef_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        namef = pyarray_to_ndarray(namef_tmp);
    }
    if (!pyarray_check(halofid_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofid = pyarray_to_ndarray(halofid_tmp);
    }
    if (!pyarray_check(centerc_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerc = pyarray_to_ndarray(centerc_tmp);
    }
    if (!pyarray_check(centerh_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        centerh = pyarray_to_ndarray(centerh_tmp);
    }
    if (!pyarray_check(oldname_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        oldname = pyarray_to_ndarray(oldname_tmp);
    }
    if (!pyarray_check(airDiamond_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        airDiamond = pyarray_to_ndarray(airDiamond_tmp);
    }
    if (!pyarray_check(n1_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        n1 = pyarray_to_ndarray(n1_tmp);
    }
    if (!pyarray_check(n2_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        n2 = pyarray_to_ndarray(n2_tmp);
    }
    if (!pyarray_check(n3_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        n3 = pyarray_to_ndarray(n3_tmp);
    }
    if (!pyarray_check(n4_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        n4 = pyarray_to_ndarray(n4_tmp);
    }
    if (!pyarray_check(normalf_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        normalf = pyarray_to_ndarray(normalf_tmp);
    }
    if (!pyarray_check(shift_tmp, NPY_DOUBLE, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        shift = pyarray_to_ndarray(shift_tmp);
    }
    if (!pyarray_check(Pbordnode_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        Pbordnode = pyarray_to_ndarray(Pbordnode_tmp);
    }
    if (!pyarray_check(Pbordface_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        Pbordface = pyarray_to_ndarray(Pbordface_tmp);
    }
    if (!pyarray_check(Px_face_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        Px_face = pyarray_to_ndarray(Px_face_tmp);
    }
    if (!pyarray_check(Py_face_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        Py_face = pyarray_to_ndarray(Py_face_tmp);
    }
    if (!pyarray_check(Pz_face_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        Pz_face = pyarray_to_ndarray(Pz_face_tmp);
    }
    if (!pyarray_check(BCdirichlet_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        BCdirichlet = pyarray_to_ndarray(BCdirichlet_tmp);
    }
    if (!pyarray_check(innerfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        innerfaces = pyarray_to_ndarray(innerfaces_tmp);
    }
    if (!pyarray_check(halofaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        halofaces = pyarray_to_ndarray(halofaces_tmp);
    }
    if (!pyarray_check(neumannfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        neumannfaces = pyarray_to_ndarray(neumannfaces_tmp);
    }
    if (!pyarray_check(dirichletfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        dirichletfaces = pyarray_to_ndarray(dirichletfaces_tmp);
    }
    if (!pyarray_check(periodicfaces_tmp, NPY_LONG, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        periodicfaces = pyarray_to_ndarray(periodicfaces_tmp);
    }
    bind_c_compute_p_gradient_3d(nd_ndim(&val_c, 0), nd_data(&val_c), nd_ndim(&v_ghost, 0), nd_data(&v_ghost), nd_ndim(&v_halo, 0), nd_data(&v_halo), nd_ndim(&v_node, 0), nd_data(&v_node), nd_ndim(&cellidf, 0), nd_ndim(&cellidf, 1), nd_data(&cellidf), nd_ndim(&nodeidf, 0), nd_ndim(&nodeidf, 1), nd_data(&nodeidf), nd_ndim(&centergf, 0), nd_ndim(&centergf, 1), nd_data(&centergf), nd_ndim(&namef, 0), nd_data(&namef), nd_ndim(&halofid, 0), nd_data(&halofid), nd_ndim(&centerc, 0), nd_ndim(&centerc, 1), nd_data(&centerc), nd_ndim(&centerh, 0), nd_ndim(&centerh, 1), nd_data(&centerh), nd_ndim(&oldname, 0), nd_data(&oldname), nd_ndim(&airDiamond, 0), nd_data(&airDiamond), nd_ndim(&n1, 0), nd_ndim(&n1, 1), nd_data(&n1), nd_ndim(&n2, 0), nd_ndim(&n2, 1), nd_data(&n2), nd_ndim(&n3, 0), nd_ndim(&n3, 1), nd_data(&n3), nd_ndim(&n4, 0), nd_ndim(&n4, 1), nd_data(&n4), nd_ndim(&normalf, 0), nd_ndim(&normalf, 1), nd_data(&normalf), nd_ndim(&shift, 0), nd_ndim(&shift, 1), nd_data(&shift), nd_ndim(&Pbordnode, 0), nd_data(&Pbordnode), nd_ndim(&Pbordface, 0), nd_data(&Pbordface), nd_ndim(&Px_face, 0), nd_data(&Px_face), nd_ndim(&Py_face, 0), nd_data(&Py_face), nd_ndim(&Pz_face, 0), nd_data(&Pz_face), nd_ndim(&BCdirichlet, 0), nd_data(&BCdirichlet), nd_ndim(&innerfaces, 0), nd_data(&innerfaces), nd_ndim(&halofaces, 0), nd_data(&halofaces), nd_ndim(&neumannfaces, 0), nd_data(&neumannfaces), nd_ndim(&dirichletfaces, 0), nd_data(&dirichletfaces), nd_ndim(&periodicfaces, 0), nd_data(&periodicfaces));
    result = Py_BuildValue("");
    free_pointer(val_c);
    free_pointer(v_ghost);
    free_pointer(v_halo);
    free_pointer(v_node);
    free_pointer(cellidf);
    free_pointer(nodeidf);
    free_pointer(centergf);
    free_pointer(namef);
    free_pointer(halofid);
    free_pointer(centerc);
    free_pointer(centerh);
    free_pointer(oldname);
    free_pointer(airDiamond);
    free_pointer(n1);
    free_pointer(n2);
    free_pointer(n3);
    free_pointer(n4);
    free_pointer(normalf);
    free_pointer(shift);
    free_pointer(Pbordnode);
    free_pointer(Pbordface);
    free_pointer(Px_face);
    free_pointer(Py_face);
    free_pointer(Pz_face);
    free_pointer(BCdirichlet);
    free_pointer(innerfaces);
    free_pointer(halofaces);
    free_pointer(neumannfaces);
    free_pointer(dirichletfaces);
    free_pointer(periodicfaces);
    return result;
}
/*........................................*/

/*........................................*/
PyObject *facetocell_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
    t_ndarray u_face = {.shape = NULL};
    t_ndarray u_c = {.shape = NULL};
    t_ndarray faceidc = {.shape = NULL};
    int64_t dim;
    PyArrayObject *u_face_tmp;
    PyArrayObject *u_c_tmp;
    PyArrayObject *faceidc_tmp;
    PyObject *dim_tmp;
    PyObject *result;
    static char *kwlist[] = {
        "u_face",
        "u_c",
        "faceidc",
        "dim",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!O", kwlist, &PyArray_Type, &u_face_tmp, &PyArray_Type, &u_c_tmp, &PyArray_Type, &faceidc_tmp, &dim_tmp))
    {
        return NULL;
    }
    if (!pyarray_check(u_face_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        u_face = pyarray_to_ndarray(u_face_tmp);
    }
    if (!pyarray_check(u_c_tmp, NPY_DOUBLE, 1, NO_ORDER_CHECK))
    {
        return NULL;
    }
    else
    {
        u_c = pyarray_to_ndarray(u_c_tmp);
    }
    if (!pyarray_check(faceidc_tmp, NPY_LONG, 2, NPY_ARRAY_C_CONTIGUOUS))
    {
        return NULL;
    }
    else
    {
        faceidc = pyarray_to_ndarray(faceidc_tmp);
    }
    if (PyIs_NativeInt(dim_tmp))
    {
        dim = PyInt64_to_Int64(dim_tmp);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "\"Argument must be native int\"");
        return NULL;
    }
    bind_c_facetocell(nd_ndim(&u_face, 0), nd_data(&u_face), nd_ndim(&u_c, 0), nd_data(&u_c), nd_ndim(&faceidc, 0), nd_ndim(&faceidc, 1), nd_data(&faceidc), dim);
    result = Py_BuildValue("");
    free_pointer(u_face);
    free_pointer(u_c);
    free_pointer(faceidc);
    return result;
}
/*........................................*/

static int exec_func(PyObject* m)
{
    return 0;
}

/*........................................*/

static PyMethodDef pyccel_functions_methods[] = {
    {
        "convert_solution",
        (PyCFunction)convert_solution_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "rhs_value_dirichlet_node",
        (PyCFunction)rhs_value_dirichlet_node_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "rhs_value_dirichlet_face",
        (PyCFunction)rhs_value_dirichlet_face_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "ghost_value_slip",
        (PyCFunction)ghost_value_slip_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "ghost_value_nonslip",
        (PyCFunction)ghost_value_nonslip_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "ghost_value_neumann",
        (PyCFunction)ghost_value_neumann_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "ghost_value_dirichlet",
        (PyCFunction)ghost_value_dirichlet_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "haloghost_value_neumann",
        (PyCFunction)haloghost_value_neumann_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "haloghost_value_dirichlet",
        (PyCFunction)haloghost_value_dirichlet_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "haloghost_value_nonslip",
        (PyCFunction)haloghost_value_nonslip_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "haloghost_value_slip",
        (PyCFunction)haloghost_value_slip_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "cell_gradient_2d",
        (PyCFunction)cell_gradient_2d_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "cell_gradient_3d",
        (PyCFunction)cell_gradient_3d_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "face_gradient_2d",
        (PyCFunction)face_gradient_2d_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "face_gradient_2d_uv",
        (PyCFunction)face_gradient_2d_uv_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "face_gradient_3d",
        (PyCFunction)face_gradient_3d_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "centertovertex_2d",
        (PyCFunction)centertovertex_2d_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "centertovertex_3d",
        (PyCFunction)centertovertex_3d_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "barthlimiter_2d",
        (PyCFunction)barthlimiter_2d_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "barthlimiter_3d",
        (PyCFunction)barthlimiter_3d_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "search_element",
        (PyCFunction)search_element_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "get_triplet_2d",
        (PyCFunction)get_triplet_2d_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "compute_2dmatrix_size",
        (PyCFunction)compute_2dmatrix_size_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "compute_3dmatrix_size",
        (PyCFunction)compute_3dmatrix_size_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "get_triplet_3d",
        (PyCFunction)get_triplet_3d_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "get_rhs_loc_2d",
        (PyCFunction)get_rhs_loc_2d_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "get_rhs_glob_2d",
        (PyCFunction)get_rhs_glob_2d_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "get_rhs_loc_3d",
        (PyCFunction)get_rhs_loc_3d_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "get_rhs_glob_3d",
        (PyCFunction)get_rhs_glob_3d_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "compute_P_gradient_2d",
        (PyCFunction)compute_P_gradient_2d_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "compute_P_gradient_3d",
        (PyCFunction)compute_P_gradient_3d_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    {
        "facetocell",
        (PyCFunction)facetocell_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    { NULL, NULL, 0, NULL}
};

/*........................................*/

static PyModuleDef_Slot pyccel_functions_slots[] = {
    {Py_mod_exec, exec_func},
    {0, NULL},
};

/*........................................*/

static struct PyModuleDef pyccel_functions_module = {
    PyModuleDef_HEAD_INIT,
    /* name of module */
    "pyccel_functions",
    /* module documentation, may be NULL */
    NULL,
    /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    0,
    pyccel_functions_methods,
    pyccel_functions_slots
};

/*........................................*/

PyMODINIT_FUNC PyInit_pyccel_functions(void)
{
    import_array();
    return PyModuleDef_Init(&pyccel_functions_module);
}

# https://sectionproperties.readthedocs.io/en/stable/

from matplotlib import pyplot as plt
from sectionproperties.analysis import Section
from sectionproperties.pre import Geometry, Material
from shapely import Polygon

def getSectionProperties(poly):
    geom = Geometry(poly)
    geom.create_mesh(mesh_sizes=50)
    sec = Section(geom)
    sec.calculate_geometric_properties()
    sec.calculate_warping_properties()
    sec.calculate_plastic_properties()
    #sec.display_results(fmt=".2f")


    A = sec.get_area() # Cross-section area

    c_x, c_y = sec.get_c()  # Elastic centroid (``cx``, ``cy``)

    q_x, q_y  = sec.get_q() # First moments of area about the global axis (``qx``, ``qy``)

    I_xx_c, I_yy_c, I_xy_c = sec.get_ic()  # Second moments of area about the centroidal axis (``ixx_c``, ``iyy_c``,``ixy_c``)

    z_xx_plus, z_xx_minus, z_yy_plus, z_yy_minus = sec.get_z() # Elastic section moduli about the centroidal axis with respect to the top and bottom fibres (``zxx_plus``, ``zxx_minus``, ``zyy_plus``, ``zyy_minus``)

    sxx, syy = sec.get_s() # Plastic section moduli about the centroidal axis (``sxx``, ``syy``)

    r_x, r_y = sec.get_rc() # Radii of gyration about the centroidal axis (``rx``, ``ry``)

    j = sec.get_j() # St. Venant torsion constant

    gamma = sec.get_gamma() # Warping constant

    properties = {"A": A,
                  "C x": c_x,
                  "C y": c_y,
                  "q x": q_x,
                  "q y": q_y,
                  "I xx": I_xx_c,
                  "I yy": I_yy_c,
                  "z xx": z_xx_plus,
                  "z yy": z_yy_plus,
                  "s xx": sxx,
                  "s yy": syy,
                  "r x": r_x,
                  "r y": r_y,
                  "J": j,
                  "Cw": gamma
                  }

    return sec, properties

# perimeter = sec.get_perimeter() # Cross-section perimeter
#
#     c_x, c_y = sec.get_c()  # Elastic centroid (``cx``, ``cy``)
#
#     q_x, q_y  = sec.get_q() # First moments of area about the global axis (``qx``, ``qy``)
#
#     I_xx_c, I_yy_c, I_xy_c = sec.get_ic()  # Second moments of area about the centroidal axis (``ixx_c``, ``iyy_c``,``ixy_c``)
#     I_11_c, I_22_c = sec.get_ip() # Second moments of area about the principal axis (``i11_c``, ``i22_c``)
#     I_xx_g, I_yy_g, I_xy_g = sec.get_ig()  # Second moments of area about the global axis (``ixx_g``, ``iyy_g``,``ixy_g``)
#
#     z_xx_plus, z_xx_minus, z_yy_plus, z_yy_minus = sec.get_z() # Elastic section moduli about the centroidal axis with respect to the top and bottom fibres (``zxx_plus``, ``zxx_minus``, ``zyy_plus``, ``zyy_minus``)
#     z_11_plus, z_11_minus, z_11_plus, z_11_minus = sec.get_zp() # Elastic section moduli about the principal axis with respect to the top and bottom fibres (``z11_plus``, ``z11_minus``, ``z22_plus``, ``z22_minus``)
#
#     sxx, syy = sec.get_s() # Plastic section moduli about the centroidal axis (``sxx``, ``syy``)
#     s11, s22 = sec.get_sp() # Plastic section moduli about the principal bending axis (``s11``, ``s22``)
#
#
#     r_x, r_y = sec.get_rc() # Radii of gyration about the centroidal axis (``rx``, ``ry``)
#     r_11, r_22 = sec.get_rp() # Radii of gyration about the principal axis (``r11``, ``r22``)
#
#     j = sec.get_j() # St. Venant torsion constant
#
#     gamma = sec.get_gamma() # Warping constant
#
#     x_se, y_se = sec.get_sc() # Centroidal axis shear centre (elasticity approach) (``x_se``, ``y_se``)
#     x11_se, y11_se = sec.get_sc_p() # Principal axis shear centre (elasticity approach) (``x11_se``, ``y22_se``)
#     x_st, y_st = sec.get_sc_t() #Centroidal axis shear centre (Trefftz's method) (``x_st``, ``y_st``)
#
#     a_sx, a_sy = sec.get_as() # Shear area for loading about the centroidal axis (``a_sx``, ``a_sy``)
#     a_s11, a_s22 = sec.get_as_p() # Shear area for loading about the princicpal bending axis (``a_s11``, ``a_s22``)
#
#     x_pc, y_pc = sec.get_pc() # Centroidal axis plastic centroid (``x_pc``, ``y_pc``)
#     x11_pc, y22_pc = sec.get_pc_p() # Principal bending axis plastic centroid (``x11_pc``, ``y22_pc``)
#
#     phi = sec.get_phi() # Principal bending axis angle








#
#
#
#
# concrete = Material(
#     name="Concrete",
#     elastic_modulus=30.1e3,  # MPa
#     poissons_ratio=0.2,
#     density=2.4e-6,          # kg/mm³
#     yield_strength=32,       # MPa
#     color="lightgrey")
#
#
# poly = Polygon([(0,0), (5,2), (3,7),(1,6)])
# geom = Geometry(geom=poly) #.plot_geometry()
#
# geom.material = concrete
#
#
# geom.create_mesh(mesh_sizes=100)
# sec = Section(geom)
# # sec.display_mesh_info()
# sec.plot_mesh()
# plt.show()
#
# sec.calculate_geometric_properties()
# sec.calculate_warping_properties()
# sec.calculate_plastic_properties()
# sec.display_results(fmt=".2f")
#
# ea = sec.get_ea()
# eixx, _, _ = sec.get_eic()
# ej = sec.get_ej()
# gj = sec.get_g_eff() / sec.get_e_eff() * ej
#
# print(f"Axial rigidity (E.A): {ea:.3e} N")
# print(f"Flexural rigidity (E.I): {eixx:.3e} N.mm2")
# print(f"Torsional rigidity (G.J): {gj:.3e} N.mm2")
#
# sec.plot_centroids()
#
# stress = sec.calculate_stress(n=10e3, mxx=10e6, vx=25e3, vy=50e3)
# stress.plot_stress(stress="vm", normalize=False, fmt="{x:.2f}")
#
# stress.plot_stress_vector(stress="vy_zxy", fmt="{x:.2f}")
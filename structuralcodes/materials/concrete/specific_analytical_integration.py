""" Code for analytical integrations of cross sections and material laws"""


# ----- Material model: ???

# ----- Material model: Parabola
# --- Cross section shape: rectangular
from typing import Iterable, List
from structuralcodes.core.base import Material


class CrossSection:
    def __init__(self):
        self.dummy_code = True
        self.height: List[float] = [200]
        self.width: List[float] = [200]

    def get_height(self) -> List[float]:
        return self.height

    def get_width(self) -> List[float]:
        return self.width


def integral_rect_xs_material_parabola(material: Material,
                                       cross_section: CrossSection,
                                       strain_top: float,
                                       strain_bot: float,
                                       alpha_d: float) -> (float, float):
    # Collecting variables from material and cross_section
    height: List[float] = cross_section.get_height()
    width: List[float] = cross_section.get_width()
    e_c2: float = 0.002  # material.eps_c2()  # assuming this value is positive
    n: float = 2.  # material.n()
    f_cd: float = 35 / 1.5 * 0.85  # material.f_cd()

    # Initializing calculation variables
    e_c_parabola_part: float
    rectangular_part: bool
    partitioned_height_vec: List[float]
    partitioned_width_vec: List[float]
    index_partition: int
    # fictional_height: float = 0.
    start_iteration: int
    force: float = 0.
    moment: float = 0.
    calculation_constant: float

    # Check if the strain in the bottom is larger than zero
    if strain_bot > 0:
        # the strain is positive -> relevant for column calculations. To be continued
        # fictional_height = 0.
        start_iteration = 0
        pass
    else:
        #fictional_height = 0.
        start_iteration = 0

    # Check if the strain values exceeds e_c2
    if strain_top > e_c2:
        # Must include the rectangular part
        e_c_parabola_part = e_c2
        rectangular_part = True
        # Need to find out where the strain == e_c2, and save this height to its own value in the List
        (height_vec, width_vec, index_partition) = segment_from_strain(
            height, width, strain_top, e_c_parabola_part, alpha_d, height[-1]
        )
    else:
        # Only the parabola part of the material model
        e_c_parabola_part = strain_top
        rectangular_part = False
        (height_vec, width_vec) = segment_from_strain_only_parabola(
            height, width, alpha_d, height[-1]
        )
        index_partition = len(height_vec) - 1


    # Integration starts
    # Parabola part goes first
    if strain_bot < e_c2:
        c: float = e_c_parabola_part / (e_c2 * (height_vec[index_partition] + 0))  # to be edited
        for i in range(start_iteration, index_partition + 1):
            x0 = 0. if i == 0 else height_vec[i - 1]
            (f_x0, mom_x0) = calculate_integral_parts_parabola_rectangular(x0, c, n)
            (f_x1, mom_x1) = calculate_integral_parts_parabola_rectangular(height_vec[i], c, n)
            width_f_cd: float = width_vec[i] * f_cd
            force += width_f_cd * (f_x1 - f_x0)
            moment += width_f_cd * (mom_x1 - mom_x0)

    # Rectangular part of the curve gets integrated
    if rectangular_part:
        height_rect: float
        force_rect: float
        leverage: float  # from alpha_d
        for i in range(index_partition + 1, len(width_vec)):
            height_rect = height_vec[i] - height_vec[i - 1]
            force_rect = f_cd * height_rect * width_vec[i]

            # summarizes the forces and moments
            force += force_rect
            leverage = height_rect / 2. + height_vec[i - 1]
            moment += force_rect * leverage

    d: float = alpha_d - moment / force  # this is the leverage in mm measured from the top
    return force, d


# --- Cross section shape: circular

# ----- Material model: Bilinear

# ----- Material model: Rectangular

# Methods
def segment_from_strain(height_iter: Iterable[float],
                        width_inp_iter: Iterable[float],
                        max_strain: float,
                        strain_at_shift: float,
                        alpha_d: float,
                        last_height: float) -> ([List[float], List[float], int]):
    height_vec: List[float] = []
    width_vec: List[float] = []
    e_0: float = 0.0
    e_end: float = max_strain
    s_height: float = last_height - alpha_d
    end_height: float = last_height

    # Linear interpolation to find the height at strain_at_shift.
    height_at_strain = (strain_at_shift - e_0) / (e_end - e_0) * (end_height - s_height)

    width_at_shift: float = 0.0
    index: int = 0
    not_already_inserted: bool = True

    for height_i, width_i in zip(height_iter, width_inp_iter):
        if height_i > height_at_strain and not_already_inserted:
            height_vec.insert(index, height_at_strain)
            width_vec.insert(index, width_i)
            not_already_inserted = False
        if height_i > s_height:
            width_vec.append(width_i)
            height_vec.append(height_i - s_height)
            index += 1
            if width_at_shift == 0.0 and height_i - s_height > height_at_strain:
                width_at_shift = width_i

    if index > 0:
        index -= 1

    if not_already_inserted:
        height_vec.insert(index, height_at_strain)
        width_vec.insert(index, width_at_shift)

    return height_vec, width_vec, index


def segment_from_strain_only_parabola(height_iter: Iterable[float],
                                    width_inp_iter: Iterable[float],
                                    alpha_d: float,
                                    last_height: float) -> (List[float], List[float]):
    height_vec: List[float] = []
    width_vec: List[float] = []
    ad_inverse: float = last_height - alpha_d

    for height_i, width_i in zip(height_iter, width_inp_iter):
        if ad_inverse <= height_i:
            height_vec.append(height_i - ad_inverse)
            width_vec.append(width_i)

    return height_vec, width_vec


def calculate_integral_parts_parabola_rectangular(x: float, c: float, n: float) -> (float, float):
    if x == 0.0 and n == 2.0:
        return 0.0, 0.0
    elif n == 2.0:
        return (
            -c * (c * x**3 / 3 - x**2),
            -c * (c * x**4 / 4 - 2 * x**3 / 3),
        )
    else:
        return (
            max(((1 - c * x) ** (n + 1) / (c * n + c)), 1e-10) + x,
            max(((1 - c * x) ** (n + 1) / (c ** 2 * (n + 1))), 1e-10) -
            max(((1 - c * x) ** (n + 2) / (c ** 2 * (n + 2))), 1e-10) +
            x ** 2 / 2,
        )


if __name__ == '__main__':
    xs: CrossSection = CrossSection()
    material: Material = Material(25000.)

    force, d = integral_rect_xs_material_parabola(material, xs, 0.0025, 0, 100.)
    print(force, d)


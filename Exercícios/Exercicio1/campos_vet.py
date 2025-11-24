from manim import *
import numpy as np


class VectorFieldsFrobenius(ThreeDScene):
    def construct(self):
        # 1. Setup 3D Environment
        axes = ThreeDAxes(
            x_range=[-3, 3, 1],
            y_range=[-3, 3, 1],
            z_range=[-1, 5, 1],
            x_length=6,
            y_length=6,
            z_length=4,
        )
        labels = axes.get_axis_labels(x_label="x", y_label="y", z_label="u")

        # 2. Define the Vector Fields
        # V1 = d/dx + x*d/du  => Vector(1, 0, x)
        def v1_func(point):
            x, y, u = point
            return np.array([1, 0, x])

        # V2 = d/dy + y*d/du => Vector(0, 1, y)
        def v2_func(point):
            x, y, u = point
            return np.array([0, 1, y])

        # 3. Create Arrow Groups
        # We manually create a grid of arrows for cleaner control in 3D
        v1_arrows = VGroup()
        v2_arrows = VGroup()

        grid_range = np.arange(-2, 2.1, 1)  # A coarse grid for clarity

        for x in grid_range:
            for y in grid_range:
                point = np.array([x, y, 0])

                # Create V1 arrow (RED)
                vect1 = v1_func(point)
                arrow1 = Arrow(start=point, end=point + 0.6 * vect1, buff=0, color=RED)
                v1_arrows.add(arrow1)

                # Create V2 arrow (BLUE)
                vect2 = v2_func(point)
                arrow2 = Arrow(start=point, end=point + 0.6 * vect2, buff=0, color=BLUE)
                v2_arrows.add(arrow2)

        # 4. Titles/Equations
        title_v1 = MathTex(r"V_1 = \partial_x + x \partial_u", color=RED).to_corner(UL)
        title_v2 = MathTex(r"V_2 = \partial_y + y \partial_u", color=BLUE).to_corner(UL)

        # 5. Animation Sequence
        self.set_camera_orientation(phi=70 * DEGREES, theta=-45 * DEGREES)
        self.add(axes, labels)
        self.wait(1)

        # Show V1
        self.add_fixed_in_frame_mobjects(title_v1)
        self.play(Write(title_v1))
        self.play(*[GrowArrow(a) for a in v1_arrows])
        self.wait(1)
        self.move_camera(phi=80 * DEGREES, theta=0 * DEGREES, run_time=2)
        self.wait(1)
        self.move_camera(phi=70 * DEGREES, theta=-45 * DEGREES, run_time=2)

        # Hide V1, Show V2
        self.add_fixed_in_frame_mobjects(title_v2)
        self.play(FadeOut(v1_arrows), TransformMatchingShapes(title_v1, title_v2))
        self.play(*[GrowArrow(a) for a in v2_arrows])
        self.wait(1)
        self.move_camera(phi=80 * DEGREES, theta=-90 * DEGREES, run_time=2)
        self.wait(1)

        # Show both
        self.move_camera(phi=60 * DEGREES, theta=-45 * DEGREES, run_time=2)
        self.play(FadeIn(v1_arrows))

        # Final rotation around the scene
        self.begin_ambient_camera_rotation(rate=0.2)
        self.wait(5)

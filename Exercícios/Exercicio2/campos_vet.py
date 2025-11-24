from manim import *
import numpy as np


class VectorFieldsExpSystemFilled(ThreeDScene):
    def construct(self):
        # 1. Setup 3D Environment
        axes = ThreeDAxes(
            x_range=[-2.5, 2.5, 1],
            y_range=[-2.5, 2.5, 1],
            z_range=[-2, 5, 1],
            x_length=5,
            y_length=5,
            z_length=6,
        )
        axes.shift(DOWN * 1.5)
        labels = axes.get_axis_labels(x_label="x", y_label="y", z_label="u")

        # 2. Define Vector Functions
        # X = d/dx + (1 + y + u) d/du
        def v1_func(point):
            x, y, z = point
            return np.array([1, 0, 1 + y + z])

        # Y = d/dy + (y + u) d/du
        def v2_func(point):
            x, y, z = point
            return np.array([0, 1, y + z])

        # 3. Create Arrow Groups
        v1_arrows = VGroup()
        v2_arrows = VGroup()

        # Grid Setup
        x_vals = np.linspace(-2, 2, 5)
        y_vals = np.linspace(-2, 2, 5)

        z_vals = np.arange(-1, 4, 1)

        for z in z_vals:
            for x in x_vals:
                for y in y_vals:
                    scene_point = axes.c2p(x, y, z)

                    vect1 = v1_func(np.array([x, y, z]))

                    vect1_norm = vect1 / np.linalg.norm(vect1)

                    arrow1 = Arrow(
                        start=scene_point,
                        end=scene_point + 0.4 * vect1_norm,
                        buff=0,
                        color=RED,
                        stroke_width=2,
                        max_tip_length_to_length_ratio=0.25,
                    )
                    v1_arrows.add(arrow1)

                    vect2 = v2_func(np.array([x, y, z]))
                    vect2_norm = vect2 / np.linalg.norm(vect2)

                    arrow2 = Arrow(
                        start=scene_point,
                        end=scene_point + 0.4 * vect2_norm,
                        buff=0,
                        color=BLUE,
                        stroke_width=2,
                        max_tip_length_to_length_ratio=0.25,
                    )
                    v2_arrows.add(arrow2)

        # 4. Text Labels
        title_v1 = MathTex(r"X = \partial_x + (1+y+u) \partial_u", color=RED).to_corner(
            UL
        )
        title_v2 = MathTex(r"Y = \partial_y + (y+u) \partial_u", color=BLUE).to_corner(
            UL
        )

        # 5. Animation
        self.set_camera_orientation(phi=75 * DEGREES, theta=-45 * DEGREES)
        self.add(axes, labels)

        # Animate X
        self.add_fixed_in_frame_mobjects(title_v1)
        self.play(Write(title_v1))
        self.play(*[GrowArrow(a) for a in v1_arrows], run_time=2, lag_ratio=0.01)
        self.wait(0.5)

        self.move_camera(phi=85 * DEGREES, theta=0 * DEGREES, run_time=2)
        self.wait(0.5)

        # Animate Y
        self.move_camera(phi=75 * DEGREES, theta=-45 * DEGREES, run_time=1.5)
        self.add_fixed_in_frame_mobjects(title_v2)
        self.play(FadeOut(v1_arrows), TransformMatchingShapes(title_v1, title_v2))
        self.play(*[GrowArrow(a) for a in v2_arrows], run_time=2, lag_ratio=0.01)

        # Show All
        self.play(FadeIn(v1_arrows))

        self.move_camera(phi=65 * DEGREES, theta=-30 * DEGREES, zoom=0.8, run_time=2)
        self.begin_ambient_camera_rotation(rate=0.2)
        self.wait(5)

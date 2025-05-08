from manim import *
from manim import config
#from manim.mobject.graphing.coordinate_systems import *
#from manim.mobject import *

class CreateCircle(Scene):
    def construct(self):
        circle = Circle()  # create a circle
        circle.set_fill(PINK, opacity=0.5)  # set the color and transparency
        self.play(Create(circle))  # show the circle on screen
        print("The first section is done, starting second section")
        self.wait(3)


        square = Square()  # create a square
        square.rotate(PI / 4)  # rotate a certain amount

        square.next_to(circle, RIGHT, buff=0.5)
        
        self.play(Create(square))  # animate the creation of the square
        self.play(Transform(square, circle))  # interpolate the square into the circle
        self.play(FadeOut(square))  # fade out animation



class SquareToCircle(Scene):
    def construct(self):
        circle = Circle()  # create a circle
        circle.set_fill(PINK, opacity=0.5)  # set color and transparency

        square = Square()  # create a square
        square.rotate(PI / 4)  # rotate a certain amount

        
        self.play(Create(square))  # animate the creation of the square
        self.play(Transform(square, circle))  # interpolate the square into the circle
        self.play(FadeOut(square))  # fade out animation

class SquareAndCircle(Scene):
    def construct(self):
        circle = Circle()  # create a circle
        circle.set_fill(PINK, opacity=0.5)  # set the color and transparency

        square = Square()  # create a square
        square.set_fill(BLUE, opacity=0.5)  # set the color and transparency

        square.next_to(circle, RIGHT, buff=0.5)  # set the position
        self.play(Create(circle), Create(square))  # show the shapes on screen

class AnimatedSquareToCircle(Scene):
    def construct(self):
        circle = Circle()  # create a circle
        square = Square()  # create a square

        self.play(Create(square))  # show the square on screen
        self.play(square.animate.rotate(PI / 4))  # rotate the square
        self.play(
            ReplacementTransform(square, circle)
        )  # transform the square into a circle
        self.play(
            circle.animate.set_fill(PINK, opacity=0.5)
        )  # color the circle on screen

class DifferentRotations(Scene):
    def construct(self):
        left_square = Square(color=BLUE, fill_opacity=0.7).shift(2 * LEFT)
        right_square = Square(color=GREEN, fill_opacity=0.7).shift(2 * RIGHT)
        self.play(
            left_square.animate.rotate(PI), Rotate(right_square, angle=PI), run_time=2
        )
        self.wait()


class testing(ThreeDScene):
    def construct(self):
        axis = ThreeDAxes()
        self.set_camera_orientation(phi = 75 * DEGREES, theta = -45 * DEGREES)
        self.play(FadeIn(axis))
        self.wait(3)
        sphere = Sphere(radius = 3)
        sphere.set_fill(RED, opacity = 0.5)
        self.play(FadeIn(sphere))
        self.wait(2)
        self.play(FadeOut(sphere))
        self.wait(1)


class metricSpace(Scene):
    def construct(self):
        text1 = Tex(r"1. Metric Spaces\\ 2. Complete Metric Spaces", font_size = 50, 
        tex_environment="flushleft")
        
        
        #\\begin{flalign*}
        #& 1. \\hspace{0.25cm} \\text{Metric Spaces} \\ \\\\
        #& 2. \\hspace{0.25cm} \\text{Complete Metric Spaces} 
        #\\end{flalign*}''', font_size=50)
        text1.to_edge(UL)
        text1.shift(UP*0.45 + LEFT*0.2)
        self.play(Write(text1, run_time=3))
        self.wait(3)
        self.play(Unwrite(text1, run_time=1.5, reverse = False))
        self.wait(2)
        text2a = Tex(r'''Formally, a matric space is an ordered pair $(M,d)$
        where M is a set and d is a metric on M i.e., a function $\quad d\colon M\times M\to
        \mathbb{R}
        \quad$ satisfying the following axioms for all points $x,y,z \in M$''', tex_environment = "flushleft")
        
        text2a.scale(0.7)
        
       # self.play(Write(text2a, run_time = 4))
        #self.wait(5)
        
        
        bullet1 = Tex('''\\begin{flalign*}
        & 1.\\hspace{0.25cm}\\text{Non-negativity: $d(x,y)\\geq 0$, and $d(x,y)=0$ }\\\\
        & \\text{if and only if $x=y$.}\\end{flalign*}''')
        bullet2 = Tex('''2.\\hspace{0.25cm} \\text{Symmetry: $d(x,y)=d(y,x)$}''')
        bullet3 = Tex('''3.\\hspace{0.25cm} \\text{Triangle inequality: $d(x,z)\\leq d(x,y)+d(y,z)$}''') 
        text2a.to_edge(UL)
        text2a.shift(UP*0.45 + LEFT*0.2)
        bullet1.next_to(text2a, DOWN * 2, aligned_edge = LEFT)
        bullet2.next_to(bullet1, DOWN * 0.55, aligned_edge = LEFT)
        bullet3.next_to(bullet2, DOWN * 0.55, aligned_edge = LEFT)
        self.play(Write(text2a, run_time = 2))
        self.wait(5)
        self.play(Write(bullet1, run_time = 2))
        self.play(Write(bullet2, run_time = 2))
        self.play(Write(bullet3, run_time = 2))
        self.wait(1)

class HelloLaTeX(Scene):
    def construct(self):
        tex = Tex(r"$LaTeX$", font_size = 144)
        self.add(tex)


class numberLine(ThreeDScene):
    def construct(self):
        axis1 = Axes(x_range = [-5, 5, 1], 
        x_length = 10,
        y_range = [-5, 5, 1],
        y_length = 10, 
        axis_config = {"numbers_with_elongated_ticks": [-5,5],
        "numbers_to_include": [-5, 0, 5],
        "font_size": 24})
        
        

        


        v1 = Vector([1,0], color = "BLUE", tip_length = 0.2)
        v2 = Vector([1,0], color = "YELLOW", tip_length = 0.2)
        v3 = Vector([1,0], color = "RED", tip_length = 0.2)




        v1.scale(axis1.get_x_unit_size())
        v2.scale(axis1.get_x_unit_size())
        v3.scale(axis1.get_x_unit_size())


        v1.put_start_and_end_on(axis1.c2p(0, 0), axis1.c2p(4, 0))
        v2.put_start_and_end_on(axis1.c2p(0, 0) + UP * 0.25, axis1.c2p(2, 0) + UP * 0.25)
        

        


        equation1a = Tex(r"$\lvert\lvert v_{2} - v_{1} \rvert\rvert = \sqrt{\langle v_{2} - v_{1}, v_{2} - v_{1} \rangle}$")
        equation1b = Tex(r"$= d = 4 - 2$")
        
        equation1a.move_to([-2,-2,0])
        equation1b.next_to(equation1a, DOWN)

        bullet1a = Tex(r"Non-negativity: $d(x,y)\geq 0$, and $d(x,x)=0$")
        self.play(bullet1a.animate.move_to([-5,2,0], aligned_edge = LEFT))
        
        
        

    

        self.play(FadeIn(axis1, run_time = 2))
        
        

        self.play(Create(v1))
        self.play(Create(v2))
        self.wait(1)

        
        self.wait(1)

        self.play(Write(equation1a))
        self.wait(1)
        self.play(v2.animate.put_start_and_end_on(axis1.c2p(4), axis1.c2p(2)))

        v3.put_start_and_end_on(axis1.c2p(0), axis1.c2p(2))
        self.wait(0.5)
        self.play(Create(v3))
        self.play(Write(equation1b))
        
        
        
        #self.play(v2.animate.move_to(axis1.c2p(2), aligned_edge = LEFT))


        self.wait(5)

        self.play(v1.animate.put_start_and_end_on([0, 0.5, 0], axis1.c2p(1) + [0,0.5,0]),
        v2.animate.put_start_and_end_on([0,1,0], axis1.c2p(1, 0) + [0,1,0]),
        v3.animate.put_start_and_end_on([0,1.5,0],axis1.c2p(1, 0) + [0,1.5,0]))
        self.play(Unwrite(equation1b, reverse = False))
        equation1b = Tex(r"If $v_{1} = v_{2}$ then $0 = 0$")
        equation1b.next_to(equation1a, DOWN)
        self.play(Write(equation1b))
        self.wait(1)

        self.play(v1.animate.put_start_and_end_on(axis1.c2p(0, 0), axis1.c2p(1, 0)))
        self.play(v2.animate.put_start_and_end_on(axis1.c2p(1, 0),axis1.c2p(0, 0)))
        self.wait(2)
        
        






class transitionTo3D(ThreeDScene):
    def construct(self):
        axis1 = NumberLine(x_range = [-5, 5, 1], 
        unit_size = 1, 
        numbers_with_elongated_ticks = [-10,10],
        numbers_to_include = [-5, 0, 5],
        font_size = 24)
        
        discreteMetricEquation = Tex(r"""A metric space that is a discrete: \newline\newline
        $d(p,q) = \begin{cases}0, \text{ if } p=q, \\
            1, \text{ otherwise}\end{cases}$""", font_size = 70, tex_environment= "flushleft")
        discreteMetricEquation.move_to(UP * 4 + LEFT * 6)

        self.move_camera(zoom = 0.5)

        self.play(FadeIn(axis1))
        
        self.play(Write(discreteMetricEquation))



        ax_xyz = ThreeDAxes(x_range = [-5, 5, 1], axis_config = { 
        "unit_size": 1}, tips = False, x_length = 10)
        #evans phone number 847-757-5042
        x_label = ax_xyz.get_x_axis_label(Tex("x"))
        y_label = ax_xyz.get_y_axis_label(Tex("y")).shift(UP * 1.8)
        print(x_label)
        
        self.move_camera(focal_distance = 10000)
        
        self.play(FadeIn(ax_xyz), FadeIn(x_label), FadeIn(y_label))

        dotG = Dot3D(point = ax_xyz.c2p(0, 0, 0), radius = 0.15)
        
        colors = [RED, GREEN, BLUE, PURPLE]

        dots = VGroup(
            *[dotG.copy().set_color(colors[i]) for i in range(4)]

        )
        
        

        lines = VGroup(
            #*[Line(dots[i].get_center(), dots[i+1].get_center()) for i in range(len(dots))]
        )



        
        x1 = float(3 ** (1/2))


        
        self.add(dots[0])
        
        

        self.wait(1)
        self.play(dots[1].animate.move_to(ax_xyz.c2p(x1, 0, 0)))
        
        lines += Line(dots[0].get_center(), dots[1].get_center())
        self.play(Create(lines[0]))

        self.wait(5)
        
        dots[2].move_to(dots[1])
        self.play(dots[2].animate.move_to(ax_xyz.c2p(x1/2, 3/2, 0)))
        


        lines += Line(dots[1].get_center(), dots[2].get_center())
        lines += Line(dots[2].get_center(), dots[0].get_center())
        self.play(Create(lines[1]))
        self.play(Create(lines[2]))

        dots[3].move_to(dots[2])
        
        self.wait(1)
        self.remove(axis1)
        self.move_camera(phi = 75 * DEGREES, theta = 30 * DEGREES, run_time = 1, zoom = 1)

        self.play(dots[3].animate.move_to(ax_xyz.c2p(x1/2, 1/2, (x1 - 7/16) ** (1/2))))
        
        lines += VGroup(
            *[Line(dots[3].get_center(), dots[i].get_center()) for i in range(3)]

        )
        
        
        
        self.play(Create(lines[3][0]), Create(lines[3][1]), Create(lines[3][2]))


        self.move_camera(zoom = 1.5)
        

        self.wait(2)

        self.begin_ambient_camera_rotation(rate = 0.15)

        self.wait(10)

        for line in lines:
            line.reverse_direction()
        self.play(Uncreate(lines), Uncreate(dots), Unwrite(discreteMetricEquation))
        
        self.wait(1)

class threefourfivesix(Scene):
    def construct(self):
        radius = 2
        points = [[radius*np.cos(alpha*DEGREES),radius*np.sin(alpha*DEGREES),0] for alpha in [0,180]]
        dots   = VGroup(
            *[Dot(point) for point in points]
        )
        self.play(Create(dots))
        
        print(len(dots), "HELLO WORLD")
        
        """for n in range(3,12):
            self.play(*[dots[i].animate.move_to([radius*np.cos(i*2*PI/(len(dots)+1)),radius*np.sin(i*2*PI/(len(dots)+1)),0]) for i in range(len(dots))])
                
            dots += Dot([radius*np.cos(len(dots)*2*PI/(len(dots)+1)),radius*np.sin(len(dots)*2*PI/(len(dots)+1)),0])  
            self.play(Create(dots[-1]))  

            lines = VGroup()
            for i in range(len(dots)):
                line = Line(dots[i].get_center(), dots[(i+1) % len(dots)].get_center())
                lines += line
            self.play(Create(lines))
            self.wait()
            
            for line in lines:
                line.reverse_direction()
            self.play(Uncreate(lines)) """

class completeMetricSpace2(ThreeDScene):
    def construct(self):
        self.move_camera(zoom=1)
        axis = ThreeDAxes(axis_config = {"unit_size": 1, 
                                         "include_tip": False, 
                                         "numbers_with_elongated_ticks": [-3, -1.5, 0, 2 ** (1/2), 1.5, 3],
                                         "stroke_width": 0.5}, 
                                         x_range=[-3, 3,0.25], 
                                         y_range= [-1, 1, 0.25])
        self.play(Create(axis))
        self.wait(1)
        xformula = 0
        metricdots = VGroup()
        for i in range(100):
            if i == 0:
                xformula = 1
            else:
                xformula = (xformula / 2) + (1/xformula)
            
            if i <=2:
                metricdots += Dot([xformula, 0, 0])
            else:
                metricdots += Dot([xformula, 0, 0], radius = 0.000005)
            
            
            

        self.play(Create(metricdots[0]))
        self.move_camera(zoom=3, frame_center=(2 ** (1/2), 0, 0))
        self.play(Create(metricdots[1]))
        self.play(Create(metricdots[2]))
        
        print(metricdots[0].radius)
        self.play(Transform(metricdots[1], Dot(metricdots[1].get_center(), radius = 0.0005)), Transform(metricdots[2], Dot(metricdots[2].get_center(), radius = 0.0005)))
        
        self.move_camera(zoom = 100000, frame_center=(2 ** (1/2), 0, 0))
        self.wait(2)
        self.play(Create(metricdots, rate_func = rate_functions.ease_in_out_quad, run_time = 5))
        self.wait(2)
        
        
class cauchySequenceGraph(ThreeDScene):
    def construct(self):
        self.move_camera(zoom = 0.8)
        a = 3
        b = 0.3
        f = 3
        
        def envelope(t):
            return a * np.exp(-b*t)
        
        def func(t):
            return a * np.exp(-b*t) * np.sin(2*PI*f*t)
        
        
        axis2D = Axes(x_range = [-0.5, 10, 1], 
                          y_range = [-4, 4, 1],
                          tips = False).add_coordinates()
        cauchyLabels = axis2D.get_axis_labels(Tex(r"n"), Tex(r"$X_{n}$"))
        xCauchyLabel, yCauchyLabel = cauchyLabels
        xCauchyLabel.shift(RIGHT * 6.5)
        yCauchyLabel.shift(UP * 1.8)
        self.play(Create(axis2D))          
        self.wait(1)
        funcplot = axis2D.plot(function = func,
                            x_range = [0, 10, 0.01],
                            color = YELLOW)             
        
        upperEnvPlot = axis2D.plot(function = envelope,
                                    x_range = [0, 10, 0.01],
                                    color = RED)
        
        lowerEnvPlot = axis2D.plot(function = lambda x: -envelope(x),
                            x_range = [0, 10, 0.01],
                            color = RED)
        
        self.play(Create(funcplot, run_time = 4))
        self.wait(2)
        self.play(Create(lowerEnvPlot, run_time = 3), Create(upperEnvPlot, run_time = 3))
        
        dampedSinText = Tex(r"$A \cdot B^{-\lambda t} \cdot \sin(\omega t + \phi)$")
        dampedSinExplanation = Tex(r"""$A$: Dampening Factor \newline
        $B^{-\lambda t}$: Exponential Decay Factor \newline
        $\lambda$: Damping Factor \newline
        $\omega$: Angular Frequency \newline
        $\phi$: Phase Shift""", tex_environment = "flushleft", font_size = 38)
        dampedSinText.move_to(axis2D.c2p(1, 4), aligned_edge = LEFT)
        dampedSinExplanation.move_to(axis2D.c2p(5.6, 3), aligned_edge = LEFT)
        self.play(Write(dampedSinText))
        self.wait(3)
        self.play(Write(dampedSinExplanation))
        self.wait(2)
        self.play(Unwrite(dampedSinExplanation, reverse = False),
                  Unwrite(dampedSinText, reverse = False),
                  Uncreate(funcplot),
                  Uncreate(upperEnvPlot),
                  Uncreate(lowerEnvPlot),
                  FadeOut(axis2D))
        self.wait(3)

class cauchySequenceFormula(ThreeDScene):
    def construct(self):
        cauchyDefiniton = Tex(r"""A sequence \(x_1, x_2, x_3, \ldots\) 
        in a metric space \((X, d)\) is called Cauchy if
        for every positive real number \(\epsilon > 0\) there is a positive integer 
        \(N\) such that for all positive integers \(m, n > N\), $d(x_m, x_n) < \epsilon$""", tex_environment = "flushleft", font_size = 35)
        cauchyDefiniton.move_to(UP * 3.2 + LEFT * 0.6)
        bulletPoints = VGroup()
        bulletPoints += Tex(r"A metric space \((X, d)\) is complete if any of the following equivalent conditions are satisfied:", 
        font_size = 30, tex_environment = "flushleft")
        bulletPoints += Tex(r"1. Every Cauchy sequence of points in \(X\) has a limit that is also in \(X\).", 
        font_size = 30)
        bulletPoints += Tex(r"2. Every Cauchy sequence in \(X\) converges in \(X\) (that is, to some point of \(X\)).", font_size = 30)
        bulletPoints += Tex(r"""3. Every decreasing sequence of non-empty closed subsets
         of \(X\), with diameters tending to 0, has a non-empty intersection: 
         if \(F_n\) is closed and non-empty, \(F_{n+1} \subseteq F_n\) for every \(n\), 
         and \(\operatorname{diam}(F_n) \to 0\), 
         then there is a point \(x \in X\) common to all sets \(F_n\).""", tex_environment = "flushleft", font_size = 30)
        
        prev = cauchyDefiniton
        for bullet in bulletPoints:
            
            bullet.next_to(prev, aligned_edge = LEFT, direction = DOWN * 2)
            prev = bullet

        self.play(Write(cauchyDefiniton))
        self.wait(2)
        self.play(Write(bulletPoints))
        self.wait(4)
        self.play(Unwrite(cauchyDefiniton), Unwrite(bulletPoints))

        #Example in R3
        self.set_camera_orientation(0.5)
        title = Tex(r"Example in $\mathbb{R}^3$")

        axes = ThreeDAxes(x_range = (-60, 60, 10),
                          x_length = 12,
                          y_range = (-60, 60, 10),
                          y_length = 12,
                          z_range = (-60, 60, 10),
                          z_length = 12)
        self.play(FadeIn(axes))
        self.wait(0.5)
        self.move_camera(phi = 75 * DEGREES, theta = 30 * DEGREES, run_time = 1.5)
        
        self.wait(2)
        spheres = VGroup()
        
            
        
        for i in range(50, 1, -1):
            spheres += Sphere(center = axes.c2p(0, 0, 0), radius = i, fill_opacity = 1 / (i))
        
        self.play(Create(spheres, run_time = 5))
        self.wait(2)
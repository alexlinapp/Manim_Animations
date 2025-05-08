from manim import *
from scipy import * 

class graph1(ThreeDScene):
    def construct(self):
        max_t = ValueTracker(0)
        
        
        
        axes1 = ThreeDAxes(x_range = [-20, 20, 2],
                     y_range = [-20, 20, 2],
                     x_length = 20,
                     y_length = 20,
                     z_range = [-20, 20, 2],
                     z_length = 20)
        
        def func(t):
            return np.array((np.sin(t), np.cos(t), 0))
        
        parametricFunction = ParametricFunction(func, t_range = [0, 2 * PI], fill_opacity = 0).set_color(RED)
        print("function is this type", type(ParametricFunction))
        print("func variable is", type(func))
        print(TAU)
        #self.move_camera(zoom = 0.75)
        self.add(axes1)
        def parabola(x):
            return x ** 2
        parabolaPlot = axes1.plot(parabola)
        circlePlot = axes1.plot_parametric_curve(func, t_range=[0, 2 * PI])
        
        def sine3D(t):
            return np.array((np.sin(t), t, np.cos(t)))
        
        def suction3D(u, v):
           return np.array((np.cos(u), np.sin(u), v))
      
        print(suction3D(TAU, 5))
        #springplot = axes1.plot_parametric_curve(sine3D, 
        #                                         t_range = [-5 * TAU, max_t.get_value() - 5 * TAU])
        suctionSurface = Surface(lambda u, v: axes1.c2p(*suction3D(u, v)), u_range=[0, TAU], v_range=[-5, 5])
        
        axes1.plot_surface
        
        self.add(axes1)
        self.play(Create(parabolaPlot, run_time = 3), Create(circlePlot, run_time = 3))
        #self.play(self.camera.frame.animate.move_to(circlePlot))
        self.move_camera(frame_center=circlePlot)
        self.play(Create(parametricFunction))
        self.wait(2)

        self.move_camera(phi = 60 * DEGREES, theta = 30 * DEGREES, zoom = 0.75, run_time = 3)
        #self.play(Create(springplot))
       
        self.play(Create(suctionSurface))
        self.wait(2)
        
        
        
        actualPlot4 = always_redraw(lambda : VGroup(axes1.plot_parametric_curve(sine3D, 
                                                 t_range = [-5 * TAU, max_t.get_value() - (5 * TAU), 0.1]), 
                                           Dot(axes1.c2p(*sine3D(max_t.get_value() - 5*TAU)),color=RED)))
        self.add(actualPlot4)
        self.wait(1)
        
        self.play(max_t.animate.set_value(50), rate_func=rate_functions.ease_in_out_quart, run_time = 10)
        print(max_t.get_value())

class graph2(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes(x_range = [-10, 10, 2],
                          y_range = [-10, 10, 2],
                          z_range = [-10, 10, 2],
                          tips = False)
        labels = axes.get_axis_labels("x", "y")
        self.play(Create(axes), Create(labels))
        self.move_camera(phi = 30 * DEGREES, theta = 0 * DEGREES, run_time = 2)
        self.wait(3)




class ParaSurface(ThreeDScene):
    def func(self, u, v):
        return np.array([np.cos(u) * np.cos(v), np.cos(u) * np.sin(v), u])

    def construct(self):
        axes = ThreeDAxes(x_range=[-4,4], x_length=8)
        surface = Surface(
            lambda u, v: axes.c2p(*self.func(u, v)),
            u_range=[-PI, PI],
            v_range=[0, TAU],
            resolution=8,
            fill_color=RED,
            checkerboard_colors=None
        )
        self.set_camera_orientation(theta=70 * DEGREES, phi=75 * DEGREES)
        self.add(axes, surface)



class exponentialDecay3D(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes(x_range=[-10, 10, 2], 
                          x_length=20,
                          y_range=[-10, 10, 2],
                          y_length=20,
                          z_range= [-10, 10, 2],
                          z_length=20)
        self.play(Create(axes))
        self.move_camera(phi = 60 * DEGREES, theta = 30 * DEGREES, zoom = 0.75, run_time = 3)


        a = 3
        b = 0.3
        f = 5
        
        max_t = ValueTracker(0)
        
        
        
        
        def decay2D(t):
            return a * np.exp(-b*t)
        decay2DPlot = axes.plot(decay2D, x_range=[-10, 10], color=RED, stroke_opacity=0.3)
        self.play(Create(decay2DPlot))
        
        def sin3D(t):
            return np.array([t,decay2D(t) * np.sin(t * 2 * PI * f), decay2D(t) * np.cos(t * 2 * PI * f)])
        
        
        def wormhole(u, v):
            return np.array([u, decay2D(u) * np.sin(v * 2 * PI * f), decay2D(u) * np.cos(v * 2 * PI * f)])
        
        self.wait(2)
        finalObject = VGroup()

        """wormholeSurface = Surface(lambda u, v: axes.c2p(*wormhole(u, v)), fill_color=RED, 
                                  u_range=[-2, 10],
                                  v_range=[0, 3 * TAU],
                                  fill_opacity=0.3,
                                  checkerboard_colors=None,
                                  resolution=64)"""

        

        for i in np.linspace(0, 360, 200):
            tempPlot = decay2DPlot.copy().rotate(i * DEGREES, axis = [1, 0, 0], about_point = ORIGIN)
            finalObject += tempPlot
        self.wait(2)
        #self.play(Create(wormholeSurface))
        self.play(Transform(decay2DPlot, finalObject, run_time= 3))
        
        self.wait(1)
        
        sin3DReDrawn = always_redraw(lambda : VGroup(axes.plot_parametric_curve(sin3D, 
                                                 t_range = [-2, max_t.get_value() - 2]), 
                                           Dot3D(axes.c2p(*sin3D(max_t.get_value() - 2)),color=BLUE, radius=0.1)))
        self.add(sin3DReDrawn)
        self.play(max_t.animate.set_value(12), rate_func=rate_functions.ease_in_out_quart, run_time = 7)

class funwithCoordinateSystem(MovingCameraScene):
        def construct(self):
            axes = Axes(x_range=[-20, 20, 5],
                        x_length = 20,
                        y_range = [-20, 20, 5],
                        y_length = 20,
                         tips = False)
            self.add(axes)
            self.wait(1)
            self.play(self.camera.frame.animate.scale(2))
            
            #functions to plot
            def f1(x):
                return x ** 4 + 3 * (x ** 2) + 5 * x + 3
            
            #scipy mins
            min = optimize.minimize_scalar(f1, method = "Brent").x
            print(type(min))



            
            f1Plot = axes.plot(f1, x_range=[-2, 2], color = BLUE)
            f1DerivativePlot = axes.plot_derivative_graph(f1Plot, x_range=[-2, 2], color = RED)

            minLines = axes.get_lines_to_point(axes.c2p(min, f1(min)), color= BLUE, stroke_width = 4)
            minLabel = axes.get_graph_label(graph = f1DerivativePlot, label=Tex(rf"$({min},{f1(min)})$"), x_val=min, dot=True, direction= UR)

            self.play(Create(f1Plot, run_time = 3), Create(f1DerivativePlot, run_time = 3))
            self.play(Create(minLines), Create(minLabel))
            
            self.play(self.camera.frame.animate.move_to(axes.c2p(min, f1(min))).scale(0.5))
            self.wait(1)
            self.play(Uncreate(f1Plot, f1DerivativePlot), Uncreate(minLines, run_time = 0.5))
            self.wait(1)

class coinParadox(MovingCameraScene):
    def construct(self):
        coin1 = Circle(radius=3, color=BLUE, arc_center= [0, 0, 0])
        coin2 = Circle(radius=3, color=GREEN, arc_center = [0, 6, 0])
        self.play(Create(coin1), Create(coin2), self.camera.frame.animate.scale(2).move_to([0, 3, 0]))
        self.wait(1)
        self.play(Uncreate(coin1), Uncreate(coin2))
        self.wait(1)
        coin1 = VGroup()
        for i in np.linspace(0, 360, 19):
            arc1 = None
            print(i)
            if ((i / 20) % 2 == 0):
                arc1 = Arc(radius = 3, arc_center=[0, 6, 0], start_angle=i * DEGREES, angle = 20 * DEGREES, color= BLUE)
            else:
                arc1 = Arc(radius = 3, arc_center=[0, 6, 0], start_angle=i * DEGREES, angle = 20 * DEGREES, color = RED)
            coin1 += arc1
        self.play(Create(coin1, run_time = 3))
        self.wait(1)

class abufaliaCoin(MovingCameraScene):
    def construct(self):
        radius = 3
        coin1 = Circle(radius=radius, color=BLUE, arc_center= [0, 0, 0])
        coin2 = Circle(radius=radius, color=GREEN, arc_center = [0, 6, 0])
        self.play(Create(coin1), Create(coin2), self.camera.frame.animate.scale(2).move_to([0, radius, 0]))
        self.wait(1)
        self.play(Uncreate(coin1), Uncreate(coin2))
        self.wait(1)
        coin1 = VGroup()
        colors = [RED, BLUE]
        parts = 18
        angle_part = 360/parts * DEGREES
        for i in range(parts):
            start_angle = angle_part * i
            arc1 = Arc(radius=radius, arc_center=[0, 6, 0], start_angle=start_angle, angle=angle_part, color= colors[i%2])
            coin1 += arc1
        coin1.rotate(270 * DEGREES)
        line = VGroup()
        for i in range(parts):
            length = angle_part * radius
            start = length * i
            line1 = Line(start=[start,0,0], end=[start+length, 0,0], color=colors[i%2])
            line += line1
        self.play(Create(coin1, run_time = 3))
        self.play(coin1.animate.shift(LEFT*10+DOWN*3))
        line.move_to(coin1.get_critical_point(DOWN), aligned_edge=LEFT)
        self.play(Transform(coin1, line, path_arc=-120*DEGREES, path_arc_centers=line.get_critical_point(LEFT), run_time = 3))
        self.wait(1)
class gaussLaw(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes(x_range = [-10, 10],
                          y_range = [-10, 10],
                          z_range = [-10, 10],
                          x_length = 20,
                          y_length = 20,
                          z_length = 20,
                          z_axis_config= {"include_ticks":True},
                          axis_config = {"include_ticks":True})
        self.play(Create(axes))
        self.move_camera(phi = 60 * DEGREES, theta = 30 * DEGREES, zoom = 0.5, run_time = 2)
        self.wait(1)
        self.move_camera(zoom = 0.5, run_time = 2)


        #objects to be created
        charge1 = Dot3D(point = axes.coords_to_point(0, 0, 0), color = YELLOW)
        gaussianSurface = Sphere(charge1.get_center(), radius = 3, fill_opacity = 0.5, checkerboard_colors = None, stroke_width = 0)

        #func = lambda pos: (pos * (1 / ((pos[0] ** 2 + pos[1] ** 2 + pos[2] ** 2 + 1e-10) ** (3/2))))
        func1 = lambda pos: (100 * axes.c2p(*pos) * (1 / ((axes.c2p(*pos)[0] ** 2 + axes.c2p(*pos)[1] ** 2 + axes.c2p(*pos)[2] ** 2 + 1e-10) ** (3/2))))
        
        v1Field = ArrowVectorField(func1, x_range=[-3, 3, 1], y_range=[-3, 3, 1],z_range=[-3, 3, 1], vector_config = {})
        v1Field.fit_to_coordinate_system(axes)
        self.play(Create(v1Field))

        self.play(Create(charge1))
        self.move_camera(frame_center = charge1.get_center(), zoom = 0.8)
        self.play(Create(gaussianSurface))
        self.wait(1)
        self.begin_ambient_camera_rotation(rate = 0.3)
        self.wait(3)
        self.stop_ambient_camera_rotation()
        self.move_camera(zoom = 4, frame_center = [3, 0, 0])
        
        v1 = Arrow(start = [3, 0, 0], end = [4, 0, 0], color = YELLOW, stroke_width = 18)
        self.play(Create(v1))
        self.wait(1)
        self.play(Uncreate(v1Field), Create(v1Field.get_vector([3, 0, 0])))
        
        
        s1 = Square(side_length = 5).move_to(axes.c2p(3.5, 0, 0))
        self.wait(2)
        

        self.wait(1)
        self.move_camera(zoom = 1, run_time = 2)
        self.play(Create(v1Field.get_vector([0, 3, 0])))

        self.wait(1)
        quad1VectorField = VGroup()
        for i in np.linspace(0, PI/2, 15):
            for j in np.linspace(0, PI/2, 15):
                quad1VectorField += v1Field.get_vector([np.sin(j) * np.cos(i) * 3, np.sin(j) * np.sin(i) * 3, np.cos(j) * 3])
        self.play(Create(quad1VectorField, lag_ratio = 0))
        quad1NormalVectors = VGroup()
        for i in range(len(quad1VectorField)):
            quad1NormalVectors += quad1VectorField[i].copy()
            quad1NormalVectors[i].scale(1 / quad1VectorField[i].get_length(), about_point = quad1NormalVectors[i].get_start()).set_color(YELLOW).set_opacity(0.5)
        self.play(Create(quad1NormalVectors, lag_ratio = 0))
        #self.play(Create(quad1VectorField, lag_ratio = 0))
        self.wait(2)

class praticewithVectorfields(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes(x_range = [-10, 10],
                          y_range = [-10, 10],
                          z_range = [-10, 10],
                          x_length = 20,
                          y_length = 20,
                          z_length = 20,
                          z_axis_config= {"include_ticks":False},
                          axis_config = {"include_ticks":False})
         
        func = lambda pos: (axes.c2p(*pos) * (1 / ((axes.c2p(*pos)[0] ** 2 + axes.c2p(*pos)[1] ** 2 + axes.c2p(*pos)[2] ** 2 + 1e-10) ** (3/2))))
        #(func([-1, 2, 3]))
        
        self.play(Create(ArrowVectorField(func, x_range=[-5, 5, 1], y_range=[-5, 5, 1],z_range=[-5,5, 1])))
        self.wait(2)

class ArrowExample(Scene):
    def construct(self):
        

        # the effect of buff
        square = Square(color=MAROON_A)
        arrow_3 = Arrow(start=LEFT, end=RIGHT, buff = 0, stroke_width= 18)
        arrow_4 = Arrow(start=LEFT, end=RIGHT).next_to(arrow_3, DOWN)
        g2 = VGroup(arrow_3, arrow_4, square)
        self.play(Create(g2))
        self.wait(2)
class twoDNThreeD(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes(x_range=[-5, 5],
                          y_range=[-5, 5],
                          z_range=[-5, 5],
                          x_length=10,
                          y_length=10,
                          z_length=10)
        self.play(Create(axes))
        
        
        self.move_camera(phi = 60 * DEGREES, theta = 30 * DEGREES)
        s1 = Square()
        #self.play(Create(s1))
        
        #self.move_camera(phi = 120 * DEGREES, theta = 30 * DEGREES)
        s1.move_to([3, 0, 0])
        self.add_fixed_in_frame_mobjects(s1)
        self.wait(1)
        self.remove(s1)
        self.wait(1)
        self.play(Write(s1))
        self.wait(2)
        self.remove(s1)
        func = lambda pos: (((pos[0] * UR + pos[1] * LEFT) - pos) / 3)
        v1Field = ArrowVectorField(func, x_range = [-5, 5, 1], y_range = [-5, 5, 1], z_range = [-5, 5, 1])
        v1Field.fit_to_coordinate_system(axes)
        self.play(Create(v1Field))
        self.wait(2)
        self.move_camera(phi = 0 * DEGREES, theta = -90 * DEGREES)
        self.wait(1)
        self.play(Uncreate(v1Field))
        self.play(Create(v1Field.get_vector([2, 0, 0])))
        self.play(Create(v1Field.get_vector([-2, 0, 0])))
        self.wait(2)
        
class Test(Scene):
    def construct(self):
        Tmax = ValueTracker(6)

        def plot_updater(mob, dt):
            ax = Axes([-Tmax.get_value(),Tmax.get_value()],[-1.5,1.5])
            f = ax.plot(lambda t: np.cos(t))
            ax.add(f)
            mob.become(ax)

        ax1 = Axes()
        self.add(ax1)
        ax1.add_updater(plot_updater)
        self.play(Tmax.animate.set_value(12),run_time=4)
        a1 = VGroup()
        a2 = VGroup()
        for i in range(10):
            a1 += Square(i)
        for i in range(len(a1)):
            a2 += a1[i]
            print(a2[i].side_length)
        self.wait(1)
        v1 = Vector([2, 3])
        self.play(Create(v1))
        self.wait(1)
        v1.scale(2)
        self.wait(1)
        v1.scale(0.5)
        self.wait(1)
        self.play(v1.animate.scale(2, about_point = ORIGIN))
        self.wait(2)

class textStuff(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes(x_range=[-5, 5],
                          y_range=[-5, 5],
                          z_range=[-5, 5],
                          x_length=10,
                          y_length=10,
                          z_length=10)
                                     
        self.play(Create(axes))
        self.move_camera(theta = 30 * DEGREES, phi = 60 * DEGREES)
        text1 = Tex(r"$x^2 + y^2 + z^2 = \rho^2$")
        text1.move_to([-3, 3, 0])
        self.add_fixed_in_frame_mobjects(text1)
        self.remove(text1)
        self.play(Write(text1))
        self.wait(2)
        dotx, doty, dotr, dotz, dotp = Dot3D(axes.c2p(1, 0, 0), color = BLUE), Dot3D(axes.c2p(0, 1, 0), color = YELLOW), Dot3D(axes.c2p(1, 1, 0), color = GREEN), Dot3D(axes.c2p(0, 0, 1), color = RED), Dot3D(axes.c2p(1, 1, 1), color = PURPLE)
        x = ValueTracker(1)
        y = ValueTracker(1)
        z = ValueTracker(1)
        print(ORIGIN)
        lx, ly, lz, lr, lp = Line(ORIGIN, dotx.get_center()), Line(ORIGIN, doty.get_center()), Line(ORIGIN, dotz.get_center()), Line(ORIGIN, dotx.get_center() + doty.get_center()), Line(dotr.get_center() + dotz.get_center())
        
        def lineUpdater(mob, dt):
            if mob == lx:
                mob.become(Line(ORIGIN, dotx.get_center()))
            elif mob == ly:
                mob.become(Line(ORIGIN, doty.get_center()))
            elif mob == lz:
                mob.become(Line(ORIGIN, dotz.get_center()))
            elif mob == lr:
                mob.become(Line(ORIGIN, dotx.get_center() + doty.get_center()))
            elif mob == lp:
                mob.become(Line(ORIGIN, dotr.get_center() + dotz.get_center()))
        
        def dotUpdater(dt):
            dotx.set_x(x.get_value())
            doty.set_y(y.get_value())
            dotz.set_z(z.get_value())
            dotr.set_x(x.get_value())
            dotr.set_y(y.get_value())
            dotp.set_x(x.get_value())
            dotp.set_y(y.get_value())
            dotp.set_z(z.get_value())
            
        def dotxUpdater(mob, dt):
            mob.set_x(x.get_value())
            
            
        
            
        self.add(dotx, doty, dotz, dotr, dotp, lx, ly, lz, lr, lp)
        self.add_updater(dotUpdater)
        dotx.add_updater(dotxUpdater)
        lx.add_updater(lineUpdater)
        ly.add_updater(lineUpdater)
        lz.add_updater(lineUpdater)
        lr.add_updater(lineUpdater)
        lp.add_updater(lineUpdater)
        self.play(x.animate.set_value(3))
        self.wait(1)
        self.play(y.animate.set_value(3))
        self.wait(1)
        self.play(z.animate.set_value(3))
        self.wait(1)
        
        self.play(z.animate.set_value(-2), x.animate.set_value(-1), y.animate.set_value(1))
        self.wait(2)
        print(dotx.get_center(), dotr.get_center())

class abufaliaTextStuff(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes(x_range=[-5, 5],y_range=[-5, 5],z_range=[-5, 5])
        self.play(Create(axes))
        self.move_camera(theta = 30 * DEGREES, phi = 60 * DEGREES)
        dotx, doty, dotr, dotz, dotp = Dot3D(axes.c2p(1, 0, 0), color = BLUE), Dot3D(axes.c2p(0, 1, 0), color = YELLOW), Dot3D(axes.c2p(1, 1, 0), color = GREEN), Dot3D(axes.c2p(0, 0, 1), color = RED), Dot3D(axes.c2p(1, 1, 1), color = PURPLE)
        x = ValueTracker(1, name="x")
        y = ValueTracker(1, name="y")
        z = ValueTracker(1, name="z")
        print(ORIGIN)
        lx, ly, lz, lr, lp = Line(ORIGIN, dotx.get_center()), Line(ORIGIN, doty.get_center()), Line(ORIGIN, dotz.get_center()), Line(ORIGIN, dotx.get_center() + doty.get_center()), Line(dotr.get_center() + dotz.get_center())

        ar = Angle(lx, lr, radius = 0.8, other_angle = False)
        ar1 = Angle(lx, lr, radius = 0.8 + 5 * SMALL_BUFF, other_angle = False)
        phiAngle = ArcBetweenPoints(start = 0.8 * lp.get_unit_vector(), end = 0.8 * lz.get_unit_vector(), angle = angle_between_vectors(lp.get_unit_vector(), lz.get_unit_vector()))
        thetaTex = MathTex(r"\theta")
        phiTex = MathTex(r"\phi")
        thetaTex.rotate(90 * DEGREES)
        thetaTex.move_to(ar1.point_from_proportion(0.5))
        phiTex.rotate(90 * DEGREES).rotate(90 * DEGREES, axis = [0, 1, 0])
        

        def create_line_udpater(dot):
            def lineupdater(mob, dt):
                mob.become(Line(ORIGIN, dot.get_center()))
            return lineupdater
        
        def create_dotUpdater(axes):
            setters = [f"set_{axis.name}" for axis in axes]
            def dotUpdater(dot, dt):
                for setter, axis in zip(setters, axes):
                    getattr(dot, setter)(axis.get_value())
            return dotUpdater
        
        def create_angleUpdater(line1, line2):
            def angleUpdater(angle, dt):
                global ar1
                if isinstance(angle, Angle):
                    if dotx.get_x() >= 0:
                        if dotr.get_y() >= 0:
                            angle.become(Angle(line1, line2, radius = 0.8, other_angle = False))
                            ar1 = Angle(line1, line2, radius = 0.8 + 5 * SMALL_BUFF, other_angle = False)
                        else:
                            angle.become(Angle(line1, line2, radius = 0.8, other_angle = True))
                            ar1 = Angle(line1, line2, radius = 0.8 + 5 * SMALL_BUFF, other_angle = True)
                    else:
                        if dotr.get_y() >= 0:
                            angle.become(Angle(line1, line2, radius = 0.8, other_angle = True))
                            ar1 = Angle(line1, line2, radius = 0.8 + 5 * SMALL_BUFF, other_angle = True)
                        else:
                            angle.become(Angle(line1, line2, radius = 0.8, other_angle = False))
                            ar1 = Angle(line1, line2, radius = 0.8 + 5 * SMALL_BUFF, other_angle = False)
                elif isinstance(angle, ArcBetweenPoints):
                    phiAngle.become(ArcBetweenPoints(start = 0.8 * lp.get_unit_vector(), end = 0.8 * lz.get_unit_vector(), angle = angle_between_vectors(lp.get_unit_vector(), lz.get_unit_vector())))

            return angleUpdater
        def create_texUpdater(line1, line2):
            def texUpdater(tex, dt):
                global ar1
                tex.move_to(ar1.point_from_proportion(0.5))
                
            return texUpdater
        

        
        dotx.add_updater(create_dotUpdater([x]))
        doty.add_updater(create_dotUpdater([y]))
        dotz.add_updater(create_dotUpdater([z]))
        dotr.add_updater(create_dotUpdater([x, y]))
        dotp.add_updater(create_dotUpdater([x, y, z]))

        for mob, dot in zip([lx, ly, lz, lr, lp], [dotx, doty, dotz, dotr, dotp]):
            mob.add_updater(create_line_udpater(dot))

        ar.add_updater(create_angleUpdater(lx, lr))
        thetaTex.add_updater(create_texUpdater(lx, lr))
        phiAngle.add_updater(create_angleUpdater(None, None))



        self.add(dotx, doty, dotz, dotr, dotp, lx, ly, lz, lr, lp, ar, thetaTex, phiTex, phiAngle)
        self.play(x.animate.set_value(3))
        self.wait(1)
        self.play(y.animate.set_value(3))
        self.wait(1)
        self.play(x.animate.set_value(-1), y.animate.set_value(2), z.animate.set_value(-3))
        self.wait(1)
        self.play(x.animate.set_value(2))
        self.play(y.animate.set_value(-1), z.animate.set_value(0))
        self.wait(2)






class numpyTesting(ThreeDScene):
    def construct(self):
        a1 = np.array([3, 4, 5])
        a2 = np.array([4, 5, 6])
        a3 = a1 + a2
        print(a1, a2, a3)
class BraceAnnotation(Scene):
    def construct(self):
        dot = Dot([-2, -1, 0])
        dot2 = Dot([2, 1, 0])
        line = Line(dot.get_center(), dot2.get_center()).set_color(ORANGE)
        b1 = Brace(line)
        b1text = b1.get_text("Horizontal distance")
        b2 = Brace(line, direction=line.copy().rotate(PI / 2).get_unit_vector())
        b2text = b2.get_tex("x-x_1")
        self.add(line, dot, dot2, b1, b2, b1text, b2text)
        l3 = line.copy().rotate(PI)
        print(l3)

class MovingAngle(Scene):
    def construct(self):
        rotation_center = LEFT

        theta_tracker = ValueTracker(110)
        line1 = Line(LEFT, RIGHT)
        line_moving = Line(LEFT, RIGHT)
        line_ref = line_moving.copy()
        line_moving.rotate(
            theta_tracker.get_value() * DEGREES, about_point=rotation_center
        )
        a = Angle(line1, line_moving, radius=0.5, other_angle=False)
        tex = MathTex(r"\theta").move_to(
            Angle(
                line1, line_moving, radius=0.5 + 3 * SMALL_BUFF, other_angle=False
            ).point_from_proportion(0.5)
        )

        self.add(line1, line_moving, a, tex)
        self.wait()

        line_moving.add_updater(
            lambda x: x.become(line_ref.copy()).rotate(
                theta_tracker.get_value() * DEGREES, about_point=rotation_center
            )
        )

        a.add_updater(
            lambda x: x.become(Angle(line1, line_moving, radius=0.5, other_angle=False))
        )
        tex.add_updater(
            lambda x: x.move_to(
                Angle(
                    line1, line_moving, radius=0.5 + 3 * SMALL_BUFF, other_angle=False
                ).point_from_proportion(0.5)
            )
        )
        
        def create_testUpdater(dot, dt):
            print("THIS CODE RUNS")

        a1 = Circle()
        
        a2 = Square(side_length = 4)
        
        self.add(a1, a2)
        a1.add_updater(create_testUpdater)
        print(a1.get_center(), "from outside print")
        self.wait(1)
        def testing(dot):
            print(dot, "HELLo WORLD")
            a2 = Square(side_length = 1, color = BLUE)
            print(a2.get_left(), "insdie function")
        testing("BRUH")
        
        print(a2.get_left(), "outside fucntion")

        

        l1 = Line()
        print(l1.get_unit_vector(), "is of type:", type(l1.get_unit_vector()))
        self.play(theta_tracker.animate.set_value(40))
        self.play(theta_tracker.animate.increment_value(140))
        self.play(tex.animate.set_color(RED), run_time=0.5)
        self.play(theta_tracker.animate.set_value(350))

class ArcBetweenPointsExample(Scene):
    def construct(self):
        circle = Circle(radius=2, stroke_color=GREY)
        dot_1 = Dot(color=GREEN).move_to([2, 0, 0]).scale(0.5)
        dot_1_text = Tex("(2,0)").scale(0.5).next_to(dot_1, RIGHT).set_color(BLUE)
        dot_2 = Dot(color=GREEN).move_to([0, 2, 0]).scale(0.5)
        dot_2_text = Tex("(0,2)").scale(0.5).next_to(dot_2, UP).set_color(BLUE)
        arc= ArcBetweenPoints(start=2 * RIGHT, end=2 * UP, stroke_color=YELLOW)
        self.add(circle, dot_1, dot_2, dot_1_text, dot_2_text)
        self.play(Create(arc))
        self.wait(1)
        arc1 = ArcBetweenPoints(start = [2, 0, 0], end = [1, 0, 0], angle = 270 * DEGREES, stroke_color = BLUE)
        # in this case arc_center is overriden I believe
        self.play(Create(arc1))
        self.wait(2)

class normalArcRotation(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes(x_range=[-5, 5],
                          y_range=[-5, 5],
                          z_range = [-5, 5],
                          x_length=10,
                          y_length=10,
                          z_length=10)
        self.add(axes)
        self.move_camera(theta = 30 * DEGREES, phi = 60 * DEGREES)
        self.wait(1)
        lx, ly, lz, lr, lp= Line([0, 0, 0], [3, 0, 0], stroke_color = BLUE),Line([0, 0, 0], [0, 3, 0], stroke_color = YELLOW),Line([0, 0, 0], [0, 0, 3], stroke_color = RED), Line([0, 0, 0], [3,3,0], stroke_color = PURPLE), Line([0, 0, 0], [3,3,3], stroke_color = GREEN)
        lgroup = VGroup(lx, ly, lz, lr, lp)
        angle1 = Angle(lx, lr, radius = 0.5)
        anglezp = angle_between_vectors(lz.get_unit_vector(), lp.get_unit_vector())
        n1 = get_unit_normal(lz.get_vector(), lp.get_vector())
        print(n1)
        arc2 = ArcBetweenPoints(lz.get_unit_vector(), lp.get_unit_vector(), radius = 1, normal_vector = [1, 0, 0])
        v1 = Arrow(start = [3, 3, 2], end = [2, 4, 2])
        


        
        
        
        
        self.play(Create(lgroup, lag_ratio = 0))
        self.wait(1)
        self.play(Create(angle1))
        self.wait(1)
        self.play(Create(arc2))
        self.wait(1)
        self.play(Create(v1))
        self.wait(2)
        self.begin_ambient_camera_rotation(rate = 0.5)
        self.wait(8)


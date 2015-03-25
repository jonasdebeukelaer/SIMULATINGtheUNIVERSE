def fly():
    # Do a pseudocolor plot of u.
    #DeleteAllPlots()
    #AddPlot('Pseudocolor', 'temp')

    #AddOperator("Clip")
    c = ClipAttributes()
    c.funcType = c.Sphere  # Plane, Sphere
    c.center = (0, 0, 0)
    c.radius = 10
    c.sphereInverse = 1
    SetOperatorOptions(c)
    DrawPlots()
 
    # Create the control points for the views.
    c0 = View3DAttributes()
    c0.viewNormal = (0, 0, 1)
    c0.focus = (0, 0, 0)
    c0.viewUp = (1, 0, 0)
    c0.viewAngle = 30
    c0.parallelScale = 60
    c0.nearPlane = 17.3205
    c0.farPlane = 81.9615
    c0.perspective = 1
 
    c1 = View3DAttributes()
    c1.viewNormal = (0, 1, 0)
    c1.focus = (0, 0, 0)
    c1.viewUp = (1, 0, 0)
    c1.viewAngle = 30
    c1.parallelScale = 30
    c1.nearPlane = 17.3205
    c1.farPlane = 81.9615
    c1.perspective = 1
 
    c2 = View3DAttributes()
    c2.viewNormal = (0, 0, -1)
    c2.focus = (0, 0, 0)
    c2.viewUp = (1, 0, 0)
    c2.viewAngle = 30
    c2.parallelScale = 15
    c2.nearPlane = 17.3205
    c2.farPlane = 81.9615
    c2.perspective = 1
 
    c3 = View3DAttributes()
    c3.viewNormal = (0, -1, 0)
    c3.focus = (0, 0, 0)
    c3.viewUp = (1, 0, 0)
    c3.viewAngle = 30
    c3.parallelScale = 20
    c3.nearPlane = 17.3205
    c3.farPlane = 81.9615
    c3.perspective = 1

    c4 = View3DAttributes()
    c4.viewNormal = (0, 0, 1)
    c4.focus = (0, 0, 0.1)
    c4.viewUp = (1, 0, 0)
    c4.viewAngle = 30
    c4.parallelScale = 15
    c4.nearPlane = 17.3205
    c4.farPlane = 81.9615
    c4.perspective = 1

    c5 = View3DAttributes()
    c5.viewNormal = (0, 1, 0)
    c5.focus = (0, 0, 0)
    c5.viewUp = (1, 0, 0)
    c5.viewAngle = 30
    c5.parallelScale = 30
    c5.nearPlane = 17.3205
    c5.farPlane = 81.9615
    c5.perspective = 1

 
    # Create a tuple of camera values and x values. The x values are weights
    # that help to determine where in [0,1] the control points occur.
    cpts = (c0, c1, c2, c3, c4, c5)
    x=[]
    for i in range(6):
        x = x + [float(i) / float(5.)]
 
    # Animate the camera. Note that we use the new built-in EvalCubicSpline
    # function which takes a t value from [0,1] a tuple of t values and a tuple
    # of control points. In this case, the control points are View3DAttributes
    # objects that we are using to animate the camera but they can be any object
    # that supports +, * operators.

    s = SaveWindowAttributes()
    s.format = s.PNG
    s.fileName = "movie"
    s.width, s.height = 1024,768
    s.screenCapture = 0
    SetSaveWindowAttributes(s)

    N = TimeSliderGetNStates()

    SetTimeSliderState(0)
    for i in range(0, 2*N):
        if i < N:
            print 'loading frame ', i
            SetTimeSliderState(i)

        t = float(i) / float(2*N - 1)
        c = EvalCubicSpline(t, x, cpts)
        c.nearPlane = -34.461
        c.farPlane = 34.461
        SetView3D(c)
        # For moviemaking...
        SaveWindow()
 
fly()
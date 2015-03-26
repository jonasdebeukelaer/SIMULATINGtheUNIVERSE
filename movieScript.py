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
    c0.parallelScale = 100
    c0.nearPlane = 17.3205
    c0.farPlane = 81.9615
    c0.perspective = 1
    
    c1 = View3DAttributes()
    c1.viewNormal = (0, 1, 0)
    c1.focus = (0, 0, 0)
    c1.viewUp = (1, 0, 0)
    c1.viewAngle = 30
    c1.parallelScale = 50
    c1.nearPlane = 17.3205
    c1.farPlane = 81.9615
    c1.perspective = 1
    
    c2 = View3DAttributes()
    c2.viewNormal = (0, 0, -1)
    c2.focus = (0, 0, 0)
    c2.viewUp = (1, 0, 0)
    c2.viewAngle = 30
    c2.parallelScale = 40
    c2.nearPlane = 17.3205
    c2.farPlane = 81.9615
    c2.perspective = 1
    
    c3 = c0
    
    # Create a tuple of camera values and x values. The x values are weights
    # that help to determine where in [0,1] the control points occur.
    cpts = (c0, c1, c2, c3, c4, c5)
    x=[]
    for i in range(4):
        x = x + [float(i) / float(3.)]

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
    
    nsteps = 900
    for i in range(N):
        SetTimeSliderState(i)
        
        t = float(i) / float(nsteps - 1)
        c = EvalCubicSpline(t, x, cpts)
        c.nearPlane = -34.461
        c.farPlane = 34.461
        SetView3D(c)
        # For moviemaking...
        SaveWindow()

fly()
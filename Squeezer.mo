within ;
model Squeezer
  inner Modelica.Mechanics.MultiBody.World world(
    enableAnimation=true,
    animateWorld=false,
    animateGravity=false)
    annotation (Placement(transformation(extent={{-92,20},{-72,40}})));
  Modelica.Mechanics.MultiBody.Parts.BodyShape k1(
    r={0.007,0,0},
    r_CM={0.00092,0,0},
    m=0.04325,
    I_33=2.194e-6,
    animation=true,
    animateSphere=true)
    annotation (Placement(transformation(extent={{100,-52},{120,-32}})));
  Modelica.Mechanics.MultiBody.Joints.Revolute beta(
    n={0,0,1},
    phi(start=-0.0617138900142764496358948458001, displayUnit="rad"),
    w(start=0),
    a(start=14222.4439199541138705911625887),
    useAxisFlange=true,
    animation=false)
    annotation (Placement(transformation(extent={{70,-84},{94,-60}})));
  Modelica.Mechanics.MultiBody.Joints.Revolute theta(
    phi(start=0, displayUnit="rad"),
    w(start=0),
    a(start=-10666.8329399655854029433719415),
    animation=false,
    n={0,0,1})
    annotation (Placement(transformation(extent={{92,-16},{112,4}})));
  Modelica.Mechanics.MultiBody.Parts.BodyShape k2(
    r={-0.028,0,0},
    r_CM={-0.0115,0,0},
    m=0.00365,
    I_33=4.410e-7,
    animation=true,
    animateSphere=true)
    annotation (Placement(transformation(extent={{52,-16},{72,4}})));
  Modelica.Mechanics.MultiBody.Parts.BodyShape k5(
    r={0.04,0,0},
    r_CM={0.02308,0.00916,0},
    m=0.07050,
    I_33=1.169e-5,
    animateSphere=false)
    annotation (Placement(transformation(extent={{-64,-8},{-44,12}})));
  Modelica.Mechanics.MultiBody.Joints.Revolute delta(
    n={0,0,1},
    phi(displayUnit="rad", start=0.487364979543842550225598953530),
    w(start=0),
    a(start=0),
    useAxisFlange=false,
    animation=false)
    annotation (Placement(transformation(extent={{-90,-40},{-70,-20}})));
  Modelica.Mechanics.MultiBody.Parts.Fixed fixed(r={-0.06934,-0.00227,0},
      animation=false)
    annotation (Placement(transformation(extent={{-138,-52},{-118,-32}})));
  Modelica.Mechanics.MultiBody.Joints.Revolute epsilon(
    n={0,0,1},
    w(start=0),
    a(start=0),
    phi(start=1.23054744454982119249735015568, displayUnit="rad"),
    animation=false)
    annotation (Placement(transformation(extent={{-92,-68},{-72,-48}})));
  Modelica.Mechanics.MultiBody.Parts.BodyShape k7(
    r={0,-4.e-2,0},
    r_CM={-449.e-5,-1228.e-5,0},
    m=0.05498,
    I_33=1.912e-5,
    animateSphere=false)
    annotation (Placement(transformation(extent={{-64,-86},{-44,-66}})));
  Modelica.Mechanics.MultiBody.Parts.Fixed O(r={0,0,0})
    annotation (Placement(transformation(extent={{38,-100},{58,-80}})));
  Modelica.Mechanics.MultiBody.Joints.Revolute phi(
    n={0,0,1},
    phi(displayUnit="rad", start=0.222668390165885884674473185609),
    w(start=0),
    a(start=0),
    animation=false)
    annotation (Placement(transformation(extent={{-40,18},{-20,38}})));
  Modelica.Mechanics.MultiBody.Parts.BodyShape k4(
    r={0,-2.e-2,0},
    r_CM={0,1421.e-5 - 2e-2,0},
    m=0.00706,
    I_33=5.667e-7,
    animateSphere=false)
    annotation (Placement(transformation(extent={{-14,-2},{6,18}})));
  Modelica.Mechanics.MultiBody.Joints.RevolutePlanarLoopConstraint revolute(n={
        0,0,1}, animation=false)
    annotation (Placement(transformation(extent={{-6,-30},{14,-10}})));
  Modelica.Mechanics.MultiBody.Joints.RevolutePlanarLoopConstraint revolute1(n=
        {0,0,1}, animation=false)
    annotation (Placement(transformation(extent={{6,-62},{26,-42}})));
  Modelica.Mechanics.MultiBody.Joints.RevolutePlanarLoopConstraint revolute2(n=
        {0,0,1}, animation=false)
    annotation (Placement(transformation(extent={{8,26},{28,46}})));
  Modelica.Mechanics.MultiBody.Parts.BodyShape K3(
    r={0,-35.e-3,0},
    r_CM={1043.e-5,-1874.e-5,0},
    m=0.02373,
    I_33=5.255e-6,
    animateSphere=false)
    annotation (Placement(transformation(extent={{-12,48},{8,68}})));
  Modelica.Mechanics.MultiBody.Joints.Revolute gamma(
    phi(displayUnit="rad", start=0.455279819163070380255912382449),
    animation=false,
    n={0,0,1})
    annotation (Placement(transformation(extent={{-40,68},{-20,88}})));
  Modelica.Mechanics.MultiBody.Parts.Fixed B(r={-0.03635,0.03273,0}, animation=
        false)
    annotation (Placement(transformation(extent={{-72,68},{-52,88}})));
  Modelica.Mechanics.MultiBody.Parts.FixedTranslation k32(r={2.e-2,-18.e-3,0},
      animation=true)
    annotation (Placement(transformation(extent={{38,64},{58,84}})));
  Modelica.Mechanics.MultiBody.Forces.Spring spring(
    c=4530,
    s_unstretched=0.07785,
    m=0,
    animation=true,
    showMass=false,
    width=0.003785,
    coilWidth=0.00003)
    annotation (Placement(transformation(extent={{82,66},{102,86}})));
  Modelica.Mechanics.MultiBody.Parts.Fixed C(r={0.014,0.072,0}, animation=false)
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={158,90})));
  Modelica.Mechanics.MultiBody.Joints.Revolute Omega(
    animation=false,
    phi(displayUnit="rad", start=-0.222668390165885884674473185609),
    n={0,0,-1})
    annotation (Placement(transformation(extent={{-34,-92},{-14,-72}})));
  Modelica.Mechanics.MultiBody.Parts.BodyShape K6(
    r={0,0.02,0},
    r_CM={0,0.00579,0},
    m=0.00706,
    I_33=5.667e-7,
    animation=true,
    animateSphere=false)
    annotation (Placement(transformation(extent={{-24,-62},{-4,-42}})));
  Modelica.Mechanics.Rotational.Sources.ConstantTorque constantTorque1(
      tau_constant=0.033)
    annotation (Placement(transformation(extent={{36,-66},{56,-46}})));
equation
  connect(beta.frame_b, k1.frame_a) annotation (Line(
      points={{94,-72},{94,-72},{94,-42},{100,-42}},
      color={95,95,95},
      thickness=0.5));
  connect(k1.frame_b, theta.frame_a) annotation (Line(
      points={{120,-42},{130,-42},{130,-26},{88,-26},{88,-6},{92,-6}},
      color={95,95,95},
      thickness=0.5));

  connect(theta.frame_b, k2.frame_a) annotation (Line(
      points={{112,-6},{114,-6},{114,6},{114,16},{26,16},{26,-6},{52,-6}},
      color={95,95,95},
      thickness=0.5));
  connect(delta.frame_b, k5.frame_a) annotation (Line(
      points={{-70,-30},{-68,-30},{-68,2},{-64,2}},
      color={95,95,95},
      thickness=0.5));
  connect(fixed.frame_b, delta.frame_a) annotation (Line(
      points={{-118,-42},{-118,-30},{-90,-30}},
      color={95,95,95},
      thickness=0.5));
  connect(epsilon.frame_a, fixed.frame_b) annotation (Line(
      points={{-92,-58},{-92,-56},{-118,-56},{-118,-42}},
      color={95,95,95},
      thickness=0.5));
  connect(epsilon.frame_b, k7.frame_a) annotation (Line(
      points={{-72,-58},{-66,-58},{-66,-76},{-64,-76}},
      color={95,95,95},
      thickness=0.5));
  connect(beta.frame_a, O.frame_b) annotation (Line(
      points={{70,-72},{62,-72},{62,-90},{58,-90}},
      color={95,95,95},
      thickness=0.5));
  connect(k5.frame_b, phi.frame_a) annotation (Line(
      points={{-44,2},{-40,2},{-40,14},{-44,14},{-44,28},{-40,28}},
      color={95,95,95},
      thickness=0.5));
  connect(phi.frame_b, k4.frame_a) annotation (Line(
      points={{-20,28},{-10,28},{-10,16},{-22,16},{-22,8},{-14,8}},
      color={95,95,95},
      thickness=0.5));
  connect(k4.frame_b, revolute.frame_a) annotation (Line(
      points={{6,8},{8,8},{8,-2},{-10,-2},{-10,-20},{-6,-20}},
      color={95,95,95},
      thickness=0.5));
  connect(k2.frame_b, revolute2.frame_b) annotation (Line(
      points={{72,-6},{76,-6},{76,-36},{20,-36},{20,20},{20,22},{36,22},{36,36},
          {28,36}},
      color={95,95,95},
      thickness=0.5));
  connect(revolute1.frame_b, revolute2.frame_b) annotation (Line(
      points={{26,-52},{28,-52},{28,-36},{20,-36},{20,22},{36,22},{36,36},{28,
          36}},
      color={95,95,95},
      thickness=0.5));
  connect(revolute.frame_b, revolute2.frame_b) annotation (Line(
      points={{14,-20},{20,-20},{20,22},{36,22},{36,36},{28,36}},
      color={95,95,95},
      thickness=0.5));
  connect(revolute2.frame_a, K3.frame_b) annotation (Line(
      points={{8,36},{4,36},{4,42},{16,42},{16,58},{8,58}},
      color={95,95,95},
      thickness=0.5));
  connect(K3.frame_a, gamma.frame_b) annotation (Line(
      points={{-12,58},{-18,58},{-18,78},{-20,78}},
      color={95,95,95},
      thickness=0.5));
  connect(B.frame_b, gamma.frame_a) annotation (Line(
      points={{-52,78},{-40,78}},
      color={95,95,95},
      thickness=0.5));
  connect(gamma.frame_b, k32.frame_a) annotation (Line(
      points={{-20,78},{10,78},{10,74},{38,74}},
      color={95,95,95},
      thickness=0.5));
  connect(spring.frame_a, k32.frame_b) annotation (Line(
      points={{82,76},{72,76},{72,74},{58,74}},
      color={95,95,95},
      thickness=0.5));
  connect(k7.frame_b, Omega.frame_a) annotation (Line(
      points={{-44,-76},{-38,-76},{-38,-82},{-34,-82}},
      color={95,95,95},
      thickness=0.5));
  connect(Omega.frame_b, K6.frame_a) annotation (Line(
      points={{-14,-82},{-8,-82},{-8,-68},{-26,-68},{-26,-58},{-26,-52},{-24,
          -52}},
      color={95,95,95},
      thickness=0.5));
  connect(revolute1.frame_a, K6.frame_b) annotation (Line(
      points={{6,-52},{2,-52},{-4,-52}},
      color={95,95,95},
      thickness=0.5));
  connect(spring.frame_b, C.frame_b) annotation (Line(
      points={{102,76},{120,76},{138,76},{138,90},{148,90}},
      color={95,95,95},
      thickness=0.5));
  connect(constantTorque1.flange, beta.axis) annotation (Line(points={{56,-56},
          {66,-56},{66,-46},{82,-46},{82,-60}}, color={0,0,0}));
  annotation (uses(Modelica(version="3.2.1")), Diagram(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{100,100}})));
end Squeezer;

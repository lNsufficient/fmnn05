within ;
model Squeezer
  inner Modelica.Mechanics.MultiBody.World world
    annotation (Placement(transformation(extent={{8,-98},{28,-78}})));
  Modelica.Mechanics.MultiBody.Parts.BodyShape k1(
    r={0.007,0,0},
    r_CM={0.00092,0,0},
    m=0.04325,
    I_33=2.194e-6)
    annotation (Placement(transformation(extent={{100,-52},{120,-32}})));
  Modelica.Mechanics.MultiBody.Joints.Revolute beta(
    n={0,0,1},
    phi(start=-0.0617138900142764496358948458001, displayUnit="rad"),
    w(start=0),
    a(start=14222.4439199541138705911625887))
    annotation (Placement(transformation(extent={{70,-80},{94,-56}})));
  Modelica.Mechanics.MultiBody.Joints.Revolute theta(
    n={0,0,1},
    phi(start=0, displayUnit="rad"),
    w(start=0),
    a(start=-10666.8329399655854029433719415))
    annotation (Placement(transformation(extent={{92,-16},{112,4}})));
  Modelica.Mechanics.MultiBody.Parts.BodyShape k2(
    r={-0.028,0,0},
    r_CM={-0.0115,0,0},
    m=0.00365,
    I_33=4.410e-7)
    annotation (Placement(transformation(extent={{52,-16},{72,4}})));
  Modelica.Mechanics.MultiBody.Parts.BodyShape k5(
    r={0.04,0,0},
    r_CM={0.02308,0.00916,0},
    m=0.07050,
    I_33=1.169e-5)
    annotation (Placement(transformation(extent={{-64,-8},{-44,12}})));
  Modelica.Mechanics.MultiBody.Joints.Revolute delta(
    n={0,0,1},
    phi(displayUnit="rad", start=0.487364979543842550225598953530),
    w(start=0),
    a(start=0),
    useAxisFlange=false)
    annotation (Placement(transformation(extent={{-90,-40},{-70,-20}})));
  Modelica.Mechanics.MultiBody.Parts.Fixed fixed(r={-0.06934,-0.00227,0})
    annotation (Placement(transformation(extent={{-138,-52},{-118,-32}})));
  Modelica.Mechanics.MultiBody.Joints.Revolute epsilon(
    n={0,0,1},
    w(start=0),
    a(start=0),
    phi(start=1.23054744454982119249735015568, displayUnit="rad"))
    annotation (Placement(transformation(extent={{-92,-68},{-72,-48}})));
  Modelica.Mechanics.MultiBody.Parts.BodyShape k7(
    r={0,-4.e-2,0},
    r_CM={-449.e-5,-1228.e-5,0},
    m=0.05498,
    I_33=1.912e-5)
    annotation (Placement(transformation(extent={{-64,-86},{-44,-66}})));
  Modelica.Mechanics.MultiBody.Parts.Fixed O(r={0,0,0})
    annotation (Placement(transformation(extent={{32,-76},{52,-56}})));
  Modelica.Mechanics.MultiBody.Joints.Revolute phi(
    n={0,0,1},
    phi(displayUnit="rad", start=0.222668390165885884674473185609),
    w(start=0),
    a(start=0))
    annotation (Placement(transformation(extent={{-40,18},{-20,38}})));
  Modelica.Mechanics.MultiBody.Parts.BodyShape k4(
    r={0,-2.e-2,0},
    r_CM={0,1421.e-5 - 2e-2,0},
    m=0.00706,
    I_33=5.667e-7)
    annotation (Placement(transformation(extent={{-14,-2},{6,18}})));
equation
  connect(beta.frame_a, world.frame_b) annotation (Line(
      points={{70,-68},{76,-68},{76,-88},{28,-88}},
      color={95,95,95},
      thickness=0.5));
  connect(beta.frame_b, k1.frame_a) annotation (Line(
      points={{94,-68},{96,-68},{96,-42},{96,-42},{100,-42}},
      color={95,95,95},
      thickness=0.5));
  connect(k1.frame_b, theta.frame_a) annotation (Line(
      points={{120,-42},{130,-42},{130,-26},{88,-26},{88,-26},{88,-6},{92,-6}},

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
      points={{70,-68},{62,-68},{62,-66},{52,-66}},
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
  annotation (uses(Modelica(version="3.2.1")), Diagram(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{100,100}})));
end Squeezer;

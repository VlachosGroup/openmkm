import cantera as ct


def define_reactor(feed_temp, feed_pressure, feed_velocity, feed_composition,
                   cat_area_by_rctr_vol, cat_composition,
                   wall_temp, wall_area_by_rctr_vol,
                   length, rctr_type, nodes=None):
    """
    Define Continous Stirred Tank Reactor with a catalyst reacting surface
    by adding appropriate components  to IdealGasReactor class defined in
    cantera. Due to the nature of underlying class definitions in cantera,
    CSTR is defined as a collection of reactor, catlayst surface, inlets,
    outlets, and the corresponding reservoirs. Plug Flow Reactor is defined
    as a concatenation of CSTRs.

    :param feed_temp: Temperature of the reactant feed
    :param feed_pressure: Pressure of the reactant feed
    :param feed_velocity: Velocity of the reactant feed
    :param feed_composition: composition of the reactant feed. Use gas phase
    :param cat_area_by_rctr_vol: Area of catalyst given in inverse length
    :param cat_composition: Composition of catalyst. Use surface phase
    :param wall_temp: External temperature outside outer wall of reactor
    :param wall_area_by_rctr_vol: Area of outer wall given in inverse length
    :param length: Length of the reactor along which feed flows
    :param rctr_type: Type of the reactor. Options are 'CSTR' or 'PFR'
    :param nodes: If reactor type is 'PFR', it is simulated as N CSTRS, where
                    N is given by nodes

    :return: IdealGasRector or ReactorNet with appropriate components
    """
    gas = feed_composition
    gas.TP = feed_temp, feed_pressure

    surf = cat_composition
    surf.TP = gas.TP

    if rctr_type == 'CSTR':
        rctr_len = length
    else:
        rctr_len = length/(nodes-1)
    width = 1.0         # Arbitrary value

    surf_area = rctr_len * width
    rctr_vol = surf_area / cat_area_by_rctr_vol

    feed_mass_flow_rate = feed_velocity * gas.density * width / \
                          cat_area_by_rctr_vol

    """
    In the python code, the reactor is input to the defintions of
    reactor components.
    Despite the weird way of specifying reactor components in python,
    they get added to the reactor object in the C++ code
    """

    upstream = ct.Reservoir(gas, name='upstream')
    downstream = ct.Reservoir(gas, name='downstream')



    if rctr_type == 'CSTR':
        reactor = ct.IdealGasReactor(gas, energy='off')
        reactor.volume = rctr_vol
        rctr_surf = ct.ReactorSurface(surf, reactor, A=surf_area)
        inlet = ct.MassFlowController(upstream, reactor,
                                      mdot=feed_mass_flow_rate)
        outlet = ct.PressureController(reactor, downstream,
                                       master=inlet, K=1e-5) #K is small number
        return reactor, upstream, gas

    elif rctr_type == 'PFR':
        reactors = []
        for i in range(nodes):
            r = ct.IdealGasReactor(gas, energy='off')
            r.volume = rctr_vol
            r_srf = ct.ReactorSurface(surf, r, A=surf_area)

            if not i:
                inlet = ct.MassFlowController(upstream, r,
                                              mdot=feed_mass_flow_rate)
            else:
                inlet = ct.MassFlowController(reactors[i-1], r,
                                              mdot=feed_mass_flow_rate)
            outlet = ct.PressureController(r, downstream, master=inlet, K=1e-5)

            reactors.append(r)

        reactor_net = ct.ReactorNet(reactors)
        return reactor_net, upstream, gas


def advance_to_stead_state(reactor_net, max_steps=100000, rtol=1e-9,
                           atol=1e-21):
    # Bring the PFR (coupled CSTRs) to steady state

    res = reactor_net.advance_to_steady_state(max_steps=max_steps,
                                              residual_threshold=rtol,
                                              atol=atol,
                                              return_residuals=True)

    for i, residual in enumerate(res):
        if residual > rtol:
            print("Required tolerance not reached. ")




def parse_tube_inp():
    with open('tube.inp') as fp:
        pass


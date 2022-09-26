struct ReferenceAtmosphere
    rho0
    p0
    dTdz
end

rho0 = 1.205
p0 = 101325.0
T0 = rho0 * R_d / p0

StandardIsothermalAtmosphere = ReferenceAtmosphere(rho0, p0, 0.0)
StandardIsentropicAtmosphere = ReferenceAtmosphere(rho0, p0, -g/cp_d)
ConstantDensityAtmospere = ReferenceAtmosphere(rho0, p0, -g/R_d)

function temperature(z, prof::ReferenceAtmosphere)
    return T0 + prof.dTdz * z
end

function density(z, prof::ReferenceAtmosphere)
    if prof.dTdz == 0.0
        return self.rho0 * np.exp(
            -z
            * self.g
            * self.gas_properties.M
            / (scipy.constants.R * 1000.0 * self.T0)
        )
    else
        alpha = (
            self.g
            * self.gas_properties.M
            / (self.dTdz * scipy.constants.R * 1000.0)
        )
        return (
            self.rho0
            * np.power(self.T0, alpha + 1.0)
            * np.power(self.temp(pos), -alpha - 1.0)
        )
    end
end


"""


    def temp(self, pos):
        p = np.array(pos)
        if len(p.shape) > 1:
            z = p[-1]
        else:
            z = p
        return self.T0 + self.dTdz * z

    def rho(self, pos):
        p = np.array(pos)
        if len(p.shape) > 1:
            z = p[-1]
        else:
            z = p


    def drho_dz(self, pos):
        if self.dTdz == 0.0:
            return (
                -self.g
                * self.gas_properties.M
                / (scipy.constants.R * 1000.0 * self.T0)
                * self.rho(pos)
            )
        else:
            alpha = (
                self.g
                * self.gas_properties.M
                / (self.dTdz * scipy.constants.R * 1000.0)
            )
            return (-alpha - 1.0) * np.power(self.temp(pos), -alpha - 2.0) * self.dTdz

    def p(self, pos):
        return (
            self.rho(pos)
            * scipy.constants.R
            * 1000.0
            / self.gas_properties.M
            * self.temp(pos)
        )

    def pot_temperature(self, pos):
        """
        Calculate the potential temperature at pos.
        """
        return self.temp(pos) * np.power(
            self.p(pos) / self.p0, -self.gas_properties.kappa()
        )

    def x_velocity(self, pos):
        return 0.0

    def y_velocity(self, pos):
        return 0.0

    def lapseRate(self):
        """
        Calculate lapse rate for using in the CNS-AMR compressible code
        """
        return self.dTdz * scipy.constants.R * 1000.0 / self.gas_properties.M



"""
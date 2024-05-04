function [Priors, Mu_vec, Sigma_vec, attractor,K,dim]=load_Yaml_to_seDS(DS_name, pkg_dir)

txtlfile = strcat(pkg_dir,'/models/', DS_name,'.yml');
YamlStruct = ReadYaml(txtlfile);

K=YamlStruct.K;
dim=YamlStruct.dim;
Priors=YamlStruct.Priors;
Mu=YamlStruct.Mu;
Sigma=YamlStruct.Sigma;
attractor=YamlStruct.attractor;
attractor=attractor(1:3)';

Mu_vec=reshape(Mu,dim,K);
Sigma_vec=reshape(Sigma,dim,dim,K);

end





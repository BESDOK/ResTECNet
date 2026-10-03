function net = createResTECNet(H, W, C, opts)
%CREATERESTECNET Build the Res-TECNet dlnetwork (Sect. 3.2 / 3.5).
%   net = CREATERESTECNET(H, W, C, ...) with options
%     'numBlocks'   number of residual blocks             (default 8)
%     'numFilters'  filters per convolution               (default 64)
%     'L'           index of the latest TEC map channel   (default 12)
%     'residual'    true = Res-TECNet, false = plain CNN  (default true)
%     'globalSkip'  include T_t + R output skip           (default true)
%     'circular'    circular longitude padding (false = zero padding,
%                   ablation row "zero instead of circular") (default true)
%     'baseChannel' input channel used for the global skip (default 0 = the
%                   latest map, channel L: persistence base; Res-TECNet-A
%                   uses C = L+17, the SH-AR forecast for t+Delta: adaptive
%                   base + learned correction, Sect. 3.2)
%     'uncertainty' add the second output filter of Sect. 3.5
%                   (Res-TECNet-U): output channel 1 = T_t + R,
%                   channel 2 = log-variance map Sigma_t   (default false)
%
%   Input H x W x C (W = 72 unique meridians, Sect. 3.1), output
%   H x W x 1 (or H x W x 2 with 'uncertainty'). Circular padding is
%   realized with formattable functionLayers; convolutions pad only in
%   latitude. MATLAB R2024a (dlnetwork / functionLayer / trainnet).

arguments
    H (1,1) double
    W (1,1) double
    C (1,1) double
    opts.numBlocks (1,1) double = 8
    opts.numFilters (1,1) double = 64
    opts.L (1,1) double = 12
    opts.residual (1,1) logical = true
    opts.globalSkip (1,1) logical = true
    opts.circular (1,1) logical = true
    opts.uncertainty (1,1) logical = false
    opts.baseChannel (1,1) double = 0
end

nf = opts.numFilters; B = opts.numBlocks;

if opts.circular
    padFcn = @(X) cat(2, X(:, end, :, :), X, X(:, 1, :, :));
    padDesc = 'circular padding in longitude';
    convPad = [1 1 0 0];
else
    padFcn = @(X) X;                      % identity; conv zero-pads both dims
    padDesc = 'no-op (zero padding in longitude)';
    convPad = [1 1 1 1];
end
circPad = @(name) functionLayer(padFcn, 'Formattable', true, 'Name', name, ...
    'Description', padDesc);
convOpts = {'Padding', convPad, 'WeightsInitializer', 'he'};

layers = [
    imageInputLayer([H W C], 'Name', 'in', 'Normalization', 'none')
    circPad('stem_pad')
    convolution2dLayer(3, nf, convOpts{:}, 'Name', 'stem_conv')
    batchNormalizationLayer('Name', 'stem_bn')
    reluLayer('Name', 'stem_relu')
    ];
lg = layerGraph(layers);
prev = 'stem_relu';

for b = 1:B
    p = sprintf('b%d_', b);
    blk = [
        circPad([p 'pad1'])
        convolution2dLayer(3, nf, convOpts{:}, 'Name', [p 'conv1'])
        batchNormalizationLayer('Name', [p 'bn1'])
        reluLayer('Name', [p 'relu1'])
        circPad([p 'pad2'])
        convolution2dLayer(3, nf, convOpts{:}, 'Name', [p 'conv2'])
        ];
    lg = addLayers(lg, blk);
    lg = connectLayers(lg, prev, [p 'pad1']);
    if opts.residual
        lg = addLayers(lg, additionLayer(2, 'Name', [p 'add']));
        lg = connectLayers(lg, [p 'conv2'], [p 'add/in1']);
        lg = connectLayers(lg, prev,        [p 'add/in2']);
        lg = addLayers(lg, reluLayer('Name', [p 'relu2']));
        lg = connectLayers(lg, [p 'add'], [p 'relu2']);
    else
        lg = addLayers(lg, reluLayer('Name', [p 'relu2']));
        lg = connectLayers(lg, [p 'conv2'], [p 'relu2']);
    end
    prev = [p 'relu2'];
end

nOut = 1 + double(opts.uncertainty);
head = [
    circPad('head_pad')
    convolution2dLayer(3, nOut, convOpts{:}, 'Name', 'head_conv')
    ];
lg = addLayers(lg, head);
lg = connectLayers(lg, prev, 'head_pad');

if opts.globalSkip
    Lch = opts.L; if opts.baseChannel > 0, Lch = opts.baseChannel; end
    selLast = functionLayer(@(X) X(:, :, Lch, :), 'Formattable', true, ...
        'Name', 'sel_Tt', 'Description', 'select the base map channel for the global skip');
    lg = addLayers(lg, selLast);
    lg = connectLayers(lg, 'in', 'sel_Tt');
    if opts.uncertainty
        % channel 1: T_t + R ; channel 2: log-variance (no skip)
        skipU = functionLayer(@(R, T) cat(3, R(:, :, 1, :) + T, R(:, :, 2, :)), ...
            'Formattable', true, 'NumInputs', 2, 'Name', 'global_skip', ...
            'Description', 'global skip on the mean channel only');
        lg = addLayers(lg, skipU);
        lg = connectLayers(lg, 'head_conv', 'global_skip/in1');
        lg = connectLayers(lg, 'sel_Tt',    'global_skip/in2');
    else
        lg = addLayers(lg, additionLayer(2, 'Name', 'global_skip'));
        lg = connectLayers(lg, 'head_conv', 'global_skip/in1');
        lg = connectLayers(lg, 'sel_Tt',    'global_skip/in2');
    end
end

net = dlnetwork(lg);
end

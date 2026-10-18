import torch

from torch.nn import init

from torch import nn

import torch.nn.functional as F

import pytorch_lightning as pl

from torchmetrics import F1Score, Accuracy, MeanMetric, AUROC

class MLP(pl.LightningModule):
    def __init__(self, hparams):
        super().__init__()
        self.save_hyperparameters(hparams)

        input_size = self.hparams.input_size
        target_size = self.hparams.num_classes
        self.best_val_metric = 0.0

        self.net = self.make_linear_model(input_size, target_size)
        self.net.apply(self.init_weights)

        print(self.net)

        self.criterion = F.cross_entropy

        self.acc_train = Accuracy(task='multiclass', num_classes=self.hparams.num_classes)
        self.acc_val = Accuracy(task='multiclass', num_classes=self.hparams.num_classes)
        self.acc_test = Accuracy(task='multiclass', num_classes=self.hparams.num_classes)

        self.auc_train = AUROC(task='multiclass', num_classes=self.hparams.num_classes)
        self.auc_val = AUROC(task='multiclass', num_classes=self.hparams.num_classes)
        self.auc_test = AUROC(task='multiclass', num_classes=self.hparams.num_classes)

        self.f1_train = F1Score(task='multiclass', num_classes=self.hparams.num_classes, average='none')
        self.f1_val = F1Score(task='multiclass', num_classes=self.hparams.num_classes, average='none')
        self.f1_test = F1Score(task='multiclass', num_classes=self.hparams.num_classes, average='none')

        self.loss = MeanMetric()

    @staticmethod
    def add_model_specific_args(parent_parser):
        parser = parent_parser.add_argument_group("Model")
        parser.add_argument('--learning_rate', type=float, default=3e-6)
        parser.add_argument('--num_layers', type=int, default=2)
        parser.add_argument('--num_nodes_per_layer', type=int, default=100)
        parser.add_argument('--init_type', type=str, default='kaiming')
        parser.add_argument('--weight_decay', type=float, default=0)
        parser.add_argument('--lr_scheduler_gamma', type=float, default=1)
        parser.add_argument('--lr_scheduler_step_size', type=float, default=10)
        parser.add_argument('--dropout_prob', type=float, default=0.05)
        return parent_parser

    def make_linear_model(self, input_size, target_size):
        modules = [nn.Linear(input_size,self.hparams.num_nodes_per_layer)]

        for _ in range(self.hparams.num_layers-1):
            modules.extend([nn.ReLU(), nn.Dropout(self.hparams.dropout_prob), nn.BatchNorm1d(self.hparams.num_nodes_per_layer), nn.Linear(self.hparams.num_nodes_per_layer,self.hparams.num_nodes_per_layer)])
        modules.extend([nn.ReLU(), nn.Linear(self.hparams.num_nodes_per_layer,target_size)])
        net = nn.Sequential(*modules)
        return net

    def init_weights(self, m, init_gain=0.02):
        if isinstance(m, nn.Linear):
            if self.hparams.init_type == 'normal':
                init.normal_(m.weight.data, 0, 0.001)
            elif self.hparams.init_type == 'xavier':
                init.xavier_normal_(m.weight.data, gain=init_gain)
            elif self.hparams.init_type == 'kaiming':
                init.kaiming_normal_(m.weight.data, a=0, mode='fan_in')
            elif self.hparams.init_type == 'orthogonal':
                init.orthogonal_(m.weight.data, gain=init_gain)
            if hasattr(m, 'bias') and m.bias is not None:
                torch.nn.init.constant_(m.bias.data, 0.0)

    def test_step(self, batch, _):
        x, y = batch
        y_hat = self.net(x)
        self.acc_test(y_hat, y)
        self.auc_test(y_hat, y)
        return y

    def test_epoch_end(self, gts) -> float:
        """
        Test epoch end
        """
        test_acc = self.acc_test.compute()
        self.log('test.acc', test_acc)

        test_auc = self.auc_test.compute()
        self.log('test.auc', test_auc)

        return test_acc

    def forward(self, x): # Prediction
        pred = self.net(x)
        return pred

    def predict_step(self, batch, _):
        x, y = batch
        y_hat = self.net(x)
        return y_hat

    def training_step(self, batch, _):
        output = self.shared_step(batch)
        self.acc_train(output['preds'], output['gt'])

        self.log("train.loss", output['loss'], on_epoch=True, on_step=False)
        self.log("train.acc", self.acc_train, on_epoch=True, on_step=False)

        self.auc_train(output['preds'], output['gt'])
        self.log("train.auc", self.auc_train, on_epoch=True, on_step=False)

        return output['loss']

    def validation_step(self, batch, _):
        output = self.shared_step(batch)
        self.loss(output['loss'])
        self.log("val.loss", output['loss'], on_epoch=True, on_step=False)

        self.acc_val(output['preds'], output['gt'])
        self.auc_val(output['preds'], output['gt'])
        self.f1_val(output['preds'], output['gt'])

    def validation_epoch_end(self, _):
        if self.trainer.sanity_checking:
          return

        val_acc = self.acc_val.compute()
        self.log("val.acc", val_acc, on_epoch=True, on_step=False)

        val_auc = self.auc_val.compute()
        self.log("val.auc", val_auc, on_epoch=True, on_step=False)

        val_f1 = self.f1_val.compute()
        self.log("val.f1.0", val_f1[0], on_epoch=True, on_step=False, metric_attribute='f1_val')
        self.log("val.f1.1", val_f1[1], on_epoch=True, on_step=False, metric_attribute='f1_val')
        self.log("val.f1.2", val_f1[2], on_epoch=True, on_step=False, metric_attribute='f1_val')

        self.best_val_metric = max(val_f1[2], self.best_val_metric)

        self.auc_val.reset()
        self.acc_val.reset()
        self.f1_val.reset()

    def shared_step(self, batch):
        x, y = batch
        y_hat = self.net(x)
        l = self.criterion(y_hat, y)
        return {'loss':l,'preds':y_hat,'gt':y}

    def configure_optimizers(self):
        optimizer = torch.optim.Adam(self.parameters(), lr=self.hparams.learning_rate, weight_decay=self.hparams.weight_decay)
        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, patience=20/self.hparams.val_check_n, min_lr=self.hparams.learning_rate*0.001)
        return optimizer

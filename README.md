ComPWA
======

ComPWA is a project to provide a flexible and modular Partial-Wave-Analysis framework for various use-cases. It was initiated and is developed by the Panda Collaboration (antiProton ANnihilation at DArmstadt) at Fair (Facility for Antiproton and Ion Research). But ComPWA will not just be used for the Panda physics program, but also for various other experiments to provide a commonly used tool which is stable, efficient and provides comparable results. At the moment there are many PWA-tools on the market, most used just for specific experiments and specific physics cases, some experiments even have multiple tools. But why write the same software again and again? E.g. the model describing physical processes should stay the same independent where and how there was a measurement of the process. Using the actual same implementation of the model does not only save a lot of time, it also ensures that two experiments are comparing the same thing. The same argument holds for optimization-routines and estimation-functions. It might even allow combined fitting of different experiments instead of taking the average of the results!
The natural modularization, following the considerations above, would be to separate into experiment specific information, physics (models, formalisms), estimation how good the model fits the data and optimization of free parameters. The first considerations on this where discussed with experts from different experiments and different technologies where discussed and tested. The result of this process is the first requirement document of the new PWA-Framework.
This sketch illustrates the modular concept (click it to visit the Documentation of ComPWA): 
[![ComPWA Modules](https://github.com/ComPWA/ComPWA/wiki/fw.png)](http://compwa.github.io/ComPWA/ "ComPWA Documentation")



[![Bitdeli Badge](https://d2weczhvl823v0.cloudfront.net/MathiasMichel/compwa/trend.png)](https://bitdeli.com/free "Bitdeli Badge")
